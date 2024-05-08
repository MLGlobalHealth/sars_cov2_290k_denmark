"""Script to recreate Figure S3"""

# Author: Jacob Curran Sebastian (jacob.curran@sund.ku.dk)

from argparse import ArgumentParser
from time import time

import numpy as np
import seaborn as sns
import pandas as pd
import polars as pl

from matplotlib import pyplot as plt
from tqdm import tqdm

RANGES = [
    [0, 6],
    [6, 12],
    [12, 18],
    [18, 24],
    [24, 30],
    [30, 36],
    [36, 42],
    [42, 48],
    [48, 53],
]

REGION_NAMES = [
    "Nordjylland",
    "Midtjylland",
    "Syddanmark",
    "Hovedstaden",
    "Sjælland",
]

set2 = sns.color_palette("Set2")


def parse_args():
    """Parse options to recreate Fig S2."""
    parser = ArgumentParser(
        description="Arguments for Figure S2",
    )
    parser.add_argument(
        "--metadata-path",
        type=str,
        help="Metadata file",
        required=True,
    )
    parser.add_argument(
        "--hamming-path",
        type=str,
        help="Hamming distance file",
        required=True,
    )
    parser.add_argument(
        "--output-folder",
        type=str,
        help="Path to output folder",
    )
    return parser.parse_args()


def get_plot_tools(plot):
    leg = plot.legend_
    handles = leg.legend_handles
    labels = [t.get_text() for t in leg.get_texts()]
    return handles, labels


def figS3(hamming_path, output_folder, total_df_):
    regionskode = total_df_.REGIONSKODE.unique()
    regionskode = np.sort(regionskode[~np.isnan(regionskode)])

    total_df_pl = pl.from_pandas(total_df_)

    palette = {
        "Different": set2[1],
        "Alpha": set2[0],
        "Delta": set2[2],
        "Omicron": set2[3],
        "Other": "tab:grey",
    }

    for ra, rang in enumerate(RANGES):
        start = time()
        if ra != 8:
            figs, axs = plt.subplots(6, 5, figsize=(15, 10))
        else:
            figs, axs = plt.subplots(5, 5, figsize=(15, 10))
        for w in tqdm(range(rang[0], rang[1])):
            hamming_df_week = pl.read_csv(hamming_path + str(w) + ".csv")

            sequences_week = [col.split()[0] for col in hamming_df_week.columns]

            hamming_df_week.columns = sequences_week
            hamming_df_week.index = sequences_week

            total_df_week = total_df_pl.filter(
                [idx in (sequences_week) for idx in total_df_pl["strain"]]
            )

            for r, region in enumerate(regionskode):
                total_df_week_region = total_df_week.filter(
                    total_df_week["REGIONSKODE"] == region
                )
                region_sequences = list(total_df_week_region["strain"])
                n_region_sequences = len(region_sequences)

                hamming_df_week_region = hamming_df_week.select(
                    region_sequences
                ).filter([idx in (region_sequences) for idx in hamming_df_week.index])

                match_df = pl.DataFrame(
                    hamming_df_week_region.columns, schema=["strain"]
                )
                week_region_strain_df = match_df.join(
                    total_df_week_region, on="strain", how="inner"
                )
                week_region_strain_df = week_region_strain_df.select(
                    ["strain", "variant"]
                )
                equal_variants = hamming_df_week_region.clone()
                week_region_strains = np.array(equal_variants.columns)  # .values

                for seq in week_region_strains:
                    variant = week_region_strain_df.filter(
                        week_region_strain_df["strain"] == seq
                    )["variant"][0]
                    equal = (
                        np.array(week_region_strain_df["variant"]) == variant
                    ).astype(int)
                    variant_comparison = np.zeros(len(equal)).astype(str)
                    variant_comparison[equal == 0] = "Different"
                    if seq == "Other":
                        variant_comparison[equal == 1] = "Other"
                    else:
                        if variant.split()[0] == "Omicron":
                            variant = "Omicron"
                        variant_comparison[equal == 1] = variant

                    equal_variants.replace(seq, pl.Series(variant_comparison))

                unique_hamming_distances = hamming_df_week_region.to_numpy()[
                    np.triu_indices(n_region_sequences, k=1)
                ]
                unique_equal_variants = equal_variants.to_numpy()[
                    np.triu_indices(n_region_sequences, k=1)
                ]
                mean_hamming = np.mean(unique_hamming_distances)

                df = pd.DataFrame(
                    {
                        "Distance": unique_hamming_distances,
                        "Variant": unique_equal_variants,
                    },
                    columns=["Distance", "Variant"],
                )

                sns.histplot(
                    df,
                    ax=axs[w % 6, r],
                    x="Distance",
                    bins=40,
                    stat="density",
                    hue="Variant",
                    multiple="stack",
                    palette=palette,
                )
                if w % 6 == 0:
                    axs[w % 6, r].set_title(REGION_NAMES[r])

                axs[w % 6, r].axvline(mean_hamming, color="tab:red", linestyle="--")
                axs[w % 6, r].set_ylim([0, 0.06])
                axs[w % 6, r].set_xlim([0, 100])
                axs[w % 6, r].set_xlabel("Distance")
                axs[w % 6, r].set_ylabel("Week " + str(w), size="large")
                axs[w % 6, r].label_outer()

                if (
                    ra == 8
                ):  # The 5th window on the final plot is week 53, which doesn't exist.
                    if w % 6 == 5:
                        axs[w % 6, r].axis("off")

        handles, labels = get_plot_tools(axs[-1, -2])
        lgd = figs.legend(
            handles, labels, title="Variant", bbox_to_anchor=(1.075, 0.92)
        )
        figs.suptitle("Week " + str(rang[0]) + "-" + str(rang[1] - 1))

        for w in tqdm(range(rang[0], rang[1])):
            for r, region in enumerate(regionskode):

                axs[w % 6, r].get_legend().remove()

        plt.tight_layout()

        figs.savefig(
            f"{output_folder}/figS3_{ra}.pdf",
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        finish = time()
        print(
            "Figure "
            + str(ra)
            + " done in "
            + str(np.round(finish - start, 1))
            + " seconds!"
        )


def main():
    args = parse_args()

    print("Reading in metadata...")
    metadata = pd.read_csv(args.metadata_path)
    metadata.rename(columns={"KOM.x": "KOMKODE"}, inplace=True)

    total_df = metadata.copy()
    total_df = total_df.drop(0)
    total_df = total_df.dropna(subset=["SampleDateTime"])

    print("Generating plots...")
    figS3(args.hamming_path, args.output_folder, total_df)


if __name__ == "___main__":
    main()
