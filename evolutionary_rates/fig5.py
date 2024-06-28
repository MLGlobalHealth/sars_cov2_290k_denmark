"""Main script for Figure 5 and Figures S8-10"""

# Author: Neil Scheidwasser (neil.clow@sund.ku.dk)

import calendar

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd

# pylint: disable=unused-import
import scienceplots

# pylint: enable=unused-import

from evolutionary_rates.model import regress
from evolutionary_rates.plot import pointplot, reg_coef_plot
from utils.config import (
    REF_LEN,
    REGION_MAPPING,
    VACC_MAPPING,
    REGION_PALETTE,
    VACC_PALETTE,
    VARIANT_PALETTE,
)
from utils.plot import set_size

plt.style.use("nature")
plt.rcParams["font.sans-serif"] = "Arial"

month_order = {v: k for (k, v) in enumerate(calendar.month_abbr[1:])}


def parse_args():
    """Parse arguments for Figure 5."""
    parser = ArgumentParser(
        description="Arguments for Figure 5",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="data/sequenced_individual_detailed_metadata.csv",
        help="Input file (.csv)",
    )
    return parser.parse_args()


def process_data(filepath):
    """
    Read demographic data with information about the evolutionary
    rate of the SaRS-CoV-2 sequence of each individual

    Parameters
    ----------
    filepath : str
        Path to input data

    Returns
    -------
    df : pandas.DataFrame
        Processed data
    """
    df = pd.read_csv(filepath)

    df.dropna(subset=["month", "age_groups", "regions"], axis=0, inplace=True)

    df["month"] = df["month"].apply(lambda x: calendar.month_abbr[1:][x - 1])

    df["date"] = pd.to_datetime(df.date)

    df["PERSON_ID"] = df["PERSON_ID"].astype(int).astype(str)

    df["major_major_variant"] = df.major_variant.apply(lambda x: x.split(" ")[0])

    # Number of substitutions = branch length times sequence length
    df["n_changes"] = df["branch_lengths"] * REF_LEN

    # Use abbreviations for regions
    df["regions2"] = df["regions"].replace(REGION_MAPPING)

    df["vacc_status"] = (
        df["fully_vacc"]
        .eq("Fully Vaccinated")
        .astype(int)
        .add(df["partial_vacc"].eq("Partially Vaccinated"))
        .astype(str)
    )
    df["vacc_status2"] = df["vacc_status"].replace(VACC_MAPPING)

    df.sort_values(by="month", key=lambda x: x.map(month_order), inplace=True)

    print(f"Shape: {df.shape}")

    return df


def fig5(df):
    width, height = set_size(130, "s")

    for col in ["n_changes", "rate"]:
        for nnz in [False, True]:
            fig, axs = plt.subplots(1, 4, figsize=(width * 4, height), sharey=True)

            axs = axs.ravel()

            pointplot(
                df,
                x="age_groups",
                y=col,
                nnz=nnz,
                palette="dark:k",
                xlabel="Age",
                rotate=True,
                ax=axs[0],
            )

            pointplot(
                df,
                x="regions2",
                y=col,
                nnz=nnz,
                palette=REGION_PALETTE,
                xlabel="Region",
                ax=axs[1],
            )

            pointplot(
                df,
                x="vacc_status2",
                y=col,
                nnz=nnz,
                palette=VACC_PALETTE,
                xlabel="Vaccination status",
                ax=axs[2],
            )

            pointplot(
                df,
                x="major_major_variant",
                y=col,
                nnz=nnz,
                palette=VARIANT_PALETTE,
                xlabel="Major variant",
                ax=axs[3],
            )

            for i in range(1, 4):
                axs[i].set_ylabel("")

            fig.tight_layout()
            plt.savefig(
                f"evolutionary_rates/img/{col}_nnz={nnz}.pdf",
                format="pdf",
                bbox_inches="tight",
                pad_inches=0.1,
            )
            plt.show()


def fig_s8_to_s10(df, params_ols_rates_nnz, params_ols_rates, params_zinb):
    region_order = {v: k for (k, v) in enumerate(sorted(df.regions2.unique())[1:])}

    groups = {
        "age_groups": {
            "ref": "[0-15]",
        },
        # "month": {
        #     "order": month_order,
        #     "rotate_ticks": True,
        #     "ref": "Jan"
        # },
        "regions2": {"order": region_order, "palette": REGION_PALETTE, "ref": "H"},
        "vacc_status2": {"palette": VACC_PALETTE, "ref": "Full"},
        "major_major_variant": {"ref": "Alpha", "palette": VARIANT_PALETTE},
    }

    data = {
        "ols_rate_nnz": params_ols_rates_nnz,
        "ols_rate": params_ols_rates,
        "zinb_n_changes": params_zinb,
    }

    for data_key, params in data.items():
        model_name, col = data_key.split("_")[:2]
        nnz = "nnz" in data_key

        width, height = set_size(120, "s")

        fig, axs = plt.subplots(
            1, len(groups), figsize=(width * len(groups), height), sharey=True
        )

        for i, (key, specs) in enumerate(groups.items()):
            reg_coef_plot(
                params,
                group=key,
                ref=specs["ref"],
                palette=specs.get("palette", None),
                dodge=True,
                ax=axs[i],
            )

            if i == 0:
                axs[i].set_ylabel("Regression coef.")
            else:
                axs[i].set_ylabel("")

            if key == "regions2":
                axs[i].set_xlabel(f"Region\n(Ref. group: {specs['ref']})")

            if key == "major_major_variant":
                axs[i].set_xlabel(f"Major variant\n(Ref. group: {specs['ref']})")

            if i > 1:
                axs[i].margins(x=0.3)

            axs[i].axhline(0, color="k", linestyle="dashed", alpha=0.3)

        fig.tight_layout()
        plt.savefig(
            f"evolutionary_rates/img/reg_{model_name}_{col}_nnz={nnz}.pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.1,
        )


def main():
    args = parse_args()

    df = process_data(args.input)

    fig5(df)

    # OLS regression
    _, params_ols_rates = regress(model_name="ols", df=df, y="rate")

    # OLS regression with only non-zero rates
    _, params_ols_rates_nnz = regress(
        model_name="ols", df=df.query("rate > 0"), y="rate"
    )

    # ZINB regression
    _, params_zinb = regress(model_name="zinb", df=df, y="n_changes", regularized=True)

    fig_s8_to_s10(df, params_ols_rates_nnz, params_ols_rates, params_zinb)


if __name__ == "__main__":
    main()
