# pylint: disable=invalid-name
"""Script to recreate Figure S1"""

# Author: Jacob Curran Sebastian (jacob.curran@sund.ku.dk)

import os
import time

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.dates import MonthLocator
from jax import random, vmap

from growth_rates.config import REGION_DICT, STEM_LINEAGES
from growth_rates.gaussian_process import (
    gpr_preprocess,
    model,
    predict,
    run_inference,
)
from growth_rates.utils import (
    get_lineage_counts,
    get_lineage_occurences,
)


set2 = sns.color_palette("Set2")


def parse_args():
    """Parse options for Figure S1."""
    parser = ArgumentParser(
        description="Arguments for Figure S1",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--metadata-path",
        type=str,
        help="Path to metadata",
    )
    parser.add_argument(
        "--output-folder",
        type=str,
        help="Path to output folder",
    )
    return parser.parse_args()


def get_lineage_counts_and_probabilities(strain_df, dict_lineages):
    region_day_counts_dict = {}
    region_occurences_dict = {}
    tr_dict = {}
    nstrains = len(dict_lineages)
    for reg, code in REGION_DICT.items():
        region_df = strain_df.loc[strain_df["REGIONSKODE"] == code]
        days_min, days_max = region_df.agg([min, max])["DaysSinceStart"]
        time_ranges = np.zeros((nstrains, 2))
        dict_lineages_region = dict_lineages.copy()
        keep_indices = []
        keep_keys = []
        for k, key in enumerate(dict_lineages_region.keys()):

            val = dict_lineages_region[key]
            lineage_df = region_df.loc[region_df["lineage"].isin(val)]
            days_min, days_max = lineage_df.agg([min, max])["DaysSinceStart"]
            if np.isnan(days_min + days_max):
                continue
            keep_keys += [key]
            keep_indices += [k]
            time_ranges[k, :] = days_min, days_max
        dict_lineages_region = {key: dict_lineages_region[key] for key in keep_keys}
        time_ranges = time_ranges[keep_indices, :]

        tlist = [time_ranges[i, :] for i in range(time_ranges.shape[0])]

        t_ranges = {t: v for t, v in zip(dict_lineages_region.keys(), tlist)}
        tr_dict[reg] = t_ranges
        data_for_counts = region_df.filter(["DaysSinceStart", "lineage"]).reset_index(
            drop=True
        )

        day_counts = get_lineage_counts(data_for_counts, dict_lineages_region, t_ranges)
        region_day_counts_dict[reg] = day_counts

        occurences = get_lineage_occurences(data_for_counts, dict_lineages)
        region_occurences_dict[reg] = occurences

    return region_day_counts_dict, region_occurences_dict, tr_dict


def get_r_boot_region(
    region_day_counts_dict,
    region_occurences_dict,
    output_folder,
):
    r_boot_region = {}

    for reg in REGION_DICT:
        df_counts = region_day_counts_dict[reg].copy()
        occurences = region_occurences_dict[reg].copy()

        # Growth rates output
        r_boot = {}

        for lineage in df_counts.columns[1:]:
            X, y = occurences[lineage].T.values
            X = X.astype(float).reshape(-1, 1)
            y = y.astype(float)
            start_time = time.time()

            print(f"Loaded data for {lineage}")

            columns = ["All", lineage]
            df = df_counts[columns].dropna()
            # Days
            X0 = df.index.to_numpy()
            X0_min, X0_max = X0[[0, -1]]
            X1 = np.arange(X0_min, X0_max + 1)

            # If statement to do Gaussian Process classification

            print("Running GP regression")
            y_all = df.eval(f"All - `{lineage}`")
            y_lineage = df[lineage]

            X0 = X0

            y1 = gpr_preprocess(y_all)
            y2 = gpr_preprocess(y_lineage)

            rng_key, rng_key_predict = random.split(random.PRNGKey(0))
            samples1 = run_inference(model, rng_key, X0.flatten(), y1.flatten())
            print("First Fitting Done in " + str(time.time() - start_time) + " seconds")
            # Prediction

            vmap_args1 = (
                random.split(rng_key_predict, samples1["kernel_var"].shape[0]),
                samples1["kernel_var"],
                samples1["kernel_length"],
                samples1["kernel_noise"],
            )

            _, predictions1 = vmap(
                lambda rng_key, var, length, noise: predict(
                    rng_key, X1, y1.flatten(), X0, var, length, noise
                )
            )(*vmap_args1)

            print(
                "First Prediction Done in " + str(time.time() - start_time) + " seconds"
            )

            samples2 = run_inference(model, rng_key, X0.flatten(), y2.flatten())

            print(
                "Second Fitting Done in " + str(time.time() - start_time) + " seconds"
            )

            vmap_args2 = (
                random.split(rng_key_predict, samples2["kernel_var"].shape[0]),
                samples2["kernel_var"],
                samples2["kernel_length"],
                samples2["kernel_noise"],
            )

            _, predictions2 = vmap(
                lambda rng_key, var, length, noise: predict(
                    rng_key, X1, y2.flatten(), X0, var, length, noise
                )
            )(*vmap_args2)

            print(
                "Second Prediction Done in "
                + str(time.time() - start_time)
                + " seconds"
            )

            diffs2 = np.diff(predictions2, axis=1)
            diffs1 = np.diff(predictions1, axis=1)
            drs = diffs2 - diffs1

            r_boot[lineage] = drs

            np.savetxt(f"{output_folder}/{reg}_{lineage}.txt", drs)
        r_boot_region[reg] = r_boot

    return r_boot_region


def figS1(r_boot_region, dict_lineages, dates, tr_dict, output_folder):
    fig, axs = plt.subplots(5, 1)
    fig.set_figheight(10)
    fig.set_figwidth(15)
    for i, reg in enumerate(REGION_DICT):

        axs[i].grid(alpha=0.5)
        axs[i].set_ylim([-0.4, 0.4])
        axs[i].set_xlim([0, 365])
        axs[i].axhline(0, color="black", linestyle="--")
        axs[i].xaxis.set_ticks(np.arange(len(dates)), dates[:], rotation=75)
        axs[i].xaxis.set_major_locator(MonthLocator(np.arange(13).astype(int) + 1))
        axs[i].label_outer()
        r_boot_plot_dict = r_boot_region[reg]
        axs[i].set_title(reg)
        for j, lin in enumerate(dict_lineages.keys()):
            if lin not in r_boot_plot_dict.keys():
                continue
            drs = r_boot_plot_dict[lin]
            dr_means = np.mean(drs, axis=0)
            t_range = tr_dict[reg][lin]
            begin, end = t_range.astype(int)

            dr_lineage = dr_means[begin:end]
            end = begin + len(dr_lineage)

            dr_percentiles = np.percentile(drs, [5.0, 95.0], axis=0)
            axs[i].plot(dates[begin:end], dr_means[begin:end], color=set2[j], label=lin)
            axs[i].fill_between(
                dates[begin:end],
                dr_percentiles[0, begin:end],
                dr_percentiles[1, begin:end],
                color=set2[j],
                alpha=0.3,
            )

    axs[-1].set_xlabel("Date")
    handles, labels = axs[3].get_legend_handles_labels()
    fig.supylabel("Growth Rate")
    fig.suptitle("Relative Growth Rates of Variants in Danish Regions")

    fig.legend(handles, labels, bbox_to_anchor=(1.12, 0.94), title="Lineage")

    plt.tight_layout()

    plt.savefig(
        f"{output_folder}/growth_rates_per_variant_and_region.pdf",
        bbox_inches="tight",
    )


def main():
    args = parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    strain_df = (
        pd.read_csv(args.metadata_pth)
        .filter(
            [
                "strain",
                "scorpio_call",
                "SampleDateTime",
                "lineage",
                "variant",
                "REGIONSKODE",
            ],
            axis=1,
        )
        .dropna()
    )

    dates = strain_df.SampleDateTime.values

    all_lineages = sorted(list(set(strain_df.lineage)))
    dict_lineages = {
        k: [l for v0 in v for l in all_lineages if l.startswith(f"{v0}.") or l == v0]
        for k, v in STEM_LINEAGES.items()
    }

    final_dict = {}
    for k0, v0 in dict_lineages.items():
        almost_all = []
        for k1, v1 in dict_lineages.items():
            if k1 != k0 and len(v0) > len(v1):
                almost_all += v1
            final_dict[k0] = list(set(v0) - set(almost_all))

    dict_lineages = final_dict

    region_day_counts_dict, region_occurences_dict, tr_dict = (
        get_lineage_counts_and_probabilities(strain_df, dict_lineages)
    )

    r_boot_region = get_r_boot_region(
        region_day_counts_dict, region_occurences_dict, args.output_folder
    )

    figS1(r_boot_region, dict_lineages, dates, tr_dict, args.output_folder)


if __name__ == "__main__":
    main()
