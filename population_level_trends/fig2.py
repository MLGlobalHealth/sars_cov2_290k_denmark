"""Script to recreate Figure 2"""

# Author: Jacob Curran Sebastian (jacob.curran@sund.ku.dk)

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# pylint: disable=redefined-outer-name, invalid-name
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygam
import seaborn as sns

# pylint: disable=unused-import
import scienceplots

# pylint: enable=unused-import

from matplotlib.dates import MonthLocator
from tqdm import tqdm

from population_level_trends.data import load_data
from utils.config import N_DAYS, DK_POP

# Palettes
set2 = sns.color_palette("Set2")
rocket = sns.color_palette("rocket")

# Plotting style
plt.style.use("nature")
FONT_SIZE = 10
LW = 2  # linewidth
TXT_W = 0.025  # textwidth
TXT_H = 0.95  # textheight
plt.rcParams["font.size"] = FONT_SIZE
plt.rcParams["legend.title_fontsize"] = FONT_SIZE
plt.rcParams["axes.labelsize"] = FONT_SIZE
plt.rcParams["legend.fontsize"] = FONT_SIZE
plt.rcParams["xtick.labelsize"] = FONT_SIZE
plt.rcParams["xtick.labelsize"] = FONT_SIZE

# Other constants
T_RANGES = np.array(
    (
        [0, 227],
        [106, N_DAYS],
        [325, N_DAYS],
        [338, N_DAYS],
        [4, 84],
        [13, 144],
        [56, 210],
    )
)  # time ranges for each variant appearing (1st day - last day)


def parse_args():
    """Parse input-related options."""
    parser = ArgumentParser(
        description="Arguments for Figure 2",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-folder",
        type=str,
        help="Input folder path",
    )
    return parser.parse_args()


def plot_gam(seq_prop):
    X1 = np.arange(N_DAYS)
    x2 = (X1 % 7).astype(str)
    X = np.array([X1, x2]).T
    y = seq_prop

    # Fit a GAM
    gam_prop = pygam.GAM(pygam.s(0, n_splines=50) + pygam.s(1, n_splines=7)).fit(
        X, y
    )  #
    plt.figure(1).set_figwidth(15)

    XX = gam_prop.generate_X_grid(term=0)
    pdep_prop, _ = gam_prop.partial_dependence(term=0, X=XX, width=0.95)
    intercept_prop = gam_prop.coef_[-1]

    # XX_week = gam_prop.generate_X_grid(term=1)
    pdep1_prop, _ = gam_prop.partial_dependence(term=1, X=XX, width=0.95)

    plt.plot(
        XX[:, gam_prop.terms[0].feature],
        pdep_prop + intercept_prop + pdep1_prop,
        label="GAM",
    )
    plt.title(repr(gam_prop.terms[0]))
    plt.grid(alpha=0.5)
    plt.xlim([0, N_DAYS])
    plt.title("GAM fit to Proportions")
    plt.show()


def plot_nucleotide_diversity(nucleotide_diversity, days):
    """Fig S2a"""
    fig = plt.figure()

    fig.set_figwidth(15)

    plt.grid(alpha=0.5)
    plt.plot(
        np.arange(len(days) - 1),
        (np.array(nucleotide_diversity[1 : N_DAYS + 1])) * (10**3),
    )
    plt.xlabel("Time (Days)")
    plt.ylabel(r"$\pi \times 10^3$")
    plt.xlim([0, N_DAYS])
    plt.title("Nucleotide Diversity of Sequences in Denmark in 2021")
    plt.show()


def plot_relative_growth_rates(gp_folder, lineage_names):
    """Part of Fig S1"""
    dr_list = []

    fig = plt.figure()
    fig.set_figwidth(15)
    for i in tqdm(range(T_RANGES.shape[0])):
        t_min, t_max = T_RANGES[i, :]

        drs = np.genfromtxt(f"{gp_folder}/{lineage_names[i]}_r_boot_numpyro.csv")
        dr_list += [drs]
        dr_means = np.mean(drs, axis=0)
        dr_percentiles = np.percentile(drs, [1.0, 99.0], axis=0)

        t_arange = np.arange(t_min, (t_max))
        plt.plot(t_arange, dr_means[t_min:t_max], color=set2[i], label=lineage_names[i])
        plt.fill_between(
            t_arange,
            dr_percentiles[0, t_min:t_max],
            dr_percentiles[1, t_min:t_max],
            alpha=0.5,
            color=set2[i],
        )

    plt.axhline(0, color="black", linestyle="--")
    plt.xlim([0, N_DAYS + 1])
    plt.ylim([-0.4, 0.5])
    plt.ylabel("Daily Growth Rate")
    plt.xlabel("Date")
    plt.title("Daily Growth Rates of lineages in Denmark in 2021")
    plt.grid(alpha=0.5)
    fig.legend(loc="center right", bbox_to_anchor=(1.025, 0.5), title="Lineage:")
    plt.show()

    return dr_list


def fig2a(root_folder, lineage_names, dk_iar, lineage_counts, ax):
    """
    Number of sequences collected each day by date of testing,
    separated by lineage, together with the number of confirmed cases
    published each day by the Statens Serum Institut (SSI).
    """
    # Case data
    ssi_data = pd.read_csv(
        f"{root_folder}/covid_SSI_epi_data/"
        "03_bekraeftede_tilfaelde_doede_indlagte_pr_dag_pr_koen.csv",
        sep=";",
        encoding="unicode_escape",
    )

    ## Bekræftede tilfælde i alt = Total confirmed cases
    ## Prøvetagningsdato = Sampling data
    case_data = ssi_data.groupby("Prøvetagningsdato", as_index=False).sum(
        "Bekræftede tilfælde i alt"
    )
    case_data.query("'2022-01-01' >= Prøvetagningsdato >= '2021-01-01'", inplace=True)

    for i, ln in enumerate(["B.1.177"] + lineage_names + ["Other"]):
        ax.plot(
            dk_iar.date,
            lineage_counts[ln].iloc[:-1],
            label=ln,
            color=(["tab:purple"] + set2[:7] + ["tab:olive"])[i],
            linewidth=LW,
        )
    ax.set_yscale("log")
    ax.set_ylabel("Sequences / Cases")
    ax.plot(
        dk_iar.date,
        case_data["Bekræftede tilfælde i alt"].values[:-1],
        color="tab:grey",
        label="Cases",
        linewidth=LW,
    )
    ax.text(TXT_W, TXT_H, "a)", transform=ax.transAxes, fontsize=12, va="top")


def fig2b(proportions, lineage_names, dk_iar, ax):
    """
    Proportion of sequences
    collected each day belonging to each major variant
    """
    ax.stackplot(
        dk_iar.date,
        (np.flip(proportions.values[:, 1:], axis=1).T)[:, :-1],
        colors=set2[6::-1] + ["tab:purple", "tab:olive"],
        linewidth=LW,
        labels=(["Other", "B.1.117"] + lineage_names)[::-1],
    )

    ax.set_ylabel("Proportion of Sequences")
    ax.text(TXT_W, TXT_H, "b)", transform=ax.transAxes, fontsize=12, va="top")
    ax.set_ylim([0, 1])
    ax.locator_params(axis="y", nbins=3)


def fig2c(seq_prop, dk_iar, ax):
    """
    Infection Ascertainment Rate (IAR) obtained via back-calculation from
    hospitalisation and mortality data and the proportion of PCR-positive
    tests taken each day for which we have a WGS
    """
    X1 = np.arange(N_DAYS)
    X2 = (X1 % 7).astype(str)
    X = np.array([X1, X2]).T
    y = seq_prop

    # Fit a GAM
    gam_prop = pygam.GAM(pygam.s(0, n_splines=50) + pygam.s(1, n_splines=7)).fit(
        X, y
    )  #

    XX = gam_prop.generate_X_grid(term=0)
    pdep_prop, _ = gam_prop.partial_dependence(term=0, X=XX, width=0.95)
    intercept_prop = gam_prop.coef_[-1]

    # XX_week = gam_prop.generate_X_grid(term=1)
    pdep1_prop, _ = gam_prop.partial_dependence(term=1, X=XX, width=0.95)

    ax.plot(
        dk_iar.date, dk_iar.iar, color="tab:blue", label="IAR (modelled)", linewidth=LW
    )
    ax.plot(
        XX[:, gam_prop.terms[0].feature],
        pdep_prop + intercept_prop + pdep1_prop,
        label="Sequencing proportion",
        color=rocket[1],
        linewidth=LW,
    )

    ax.fill_between(
        dk_iar.date, dk_iar["iar_0.05"], dk_iar["iar_0.95"], alpha=0.4, color="tab:blue"
    )
    ax.text(TXT_W, TXT_H, "c)", transform=ax.transAxes, fontsize=12, va="top")
    ax.set_ylabel("IAR / Proportion")
    ax.set_ylim([0, 1])


def fig2d(root_folder, dk_iar, ax):
    """
    Proportion of the Danish population that have
    received a first and second vaccine dose over time
    """
    # First vaccine data
    first_vac = np.genfromtxt(f"{root_folder}/pop_data/first_vac_cumulative.csv")
    # Second vaccine data
    second_vac = np.genfromtxt(f"{root_folder}/pop_data/second_vac_cumulative.csv")
    first_second_vac_date_idx = 20

    ax.plot(
        dk_iar.date,
        first_vac / DK_POP,
        color=rocket[3],
        label="First Vaccination",
        linewidth=LW,
    )
    ax.plot(
        dk_iar.date.iloc[first_second_vac_date_idx:],
        second_vac / DK_POP,
        color=rocket[2],
        label="Second Vaccination",
        linewidth=LW,
    )
    ax.text(TXT_W, TXT_H, "d)", transform=ax.transAxes, fontsize=12, va="top")
    ax.set_ylim([0, 1])
    ax.set_ylabel("Proportion Vaccinated")


def fig2e(dk_iar, nucleotide_diversity, ax):
    """
    Nucleotide Diversity calculated
    for all sequences for each day
    """
    ax.text(TXT_W, TXT_H, "e)", transform=ax.transAxes, fontsize=12, va="top")

    ax.plot(
        dk_iar.date, np.array(nucleotide_diversity[1:366]), color="black", linewidth=LW
    )
    ax.set_ylabel("Nucleotide Diversity")
    ax.set_ylim([0, 0.0018])


def fig2f(lineage_names, dk_iar, dr_list, ax):
    """
    Daily relative growth rate calculated
    for each major lineage
    """
    # Get relative growth rates for each variant
    for i in range(T_RANGES.shape[0]):

        t_min, t_max = T_RANGES[i, :]

        drs = dr_list[i]

        dr_percentiles = np.percentile(drs, [1.0, 99.0], axis=0)

        ax.plot(color=set2[i], label=lineage_names[i], linewidth=LW)
        ax.fill_between(
            dk_iar.date.iloc[t_min:t_max],
            dr_percentiles[0, t_min:t_max],
            dr_percentiles[1, t_min:t_max],
            alpha=0.5,
            color=set2[i],
        )

    ax.axhline(0, linestyle="--", color="black")
    ax.set_ylim([-0.4, 0.4])
    ax.set_ylabel("Growth rate")
    ax.text(TXT_W, TXT_H, "f)", transform=ax.transAxes, fontsize=12, va="top")

    ax.set_xlabel("Date")


def fig2(
    root_folder,
    proportions,
    lineage_names,
    lineage_counts,
    dk_iar,
    seq_prop,
    nucleotide_diversity,
    dr_list,
):
    """Make Figure 2

    Parameters
    ----------
    root_folder : str
        Folder containing data of interest
    proportions : pandas.DataFrame
        Variant proportions
    lineage_names : List
        Lineage names
    lineage_counts : pandas.DataFrame
        Daily counts of each lineage
    dk_iar : pandas.DataFrame
        Denmark infection ascertainment rates (IARs)
    seq_prop : numpy.ndarray
        Proportion of PCR positives
    nucleotide_diversity : numpy.ndarray
        Nucleotide diversity
    dr_list : list
        Daily growth rates
    """
    plt_rows = 6
    fig, axs = plt.subplots(nrows=plt_rows, ncols=1, figsize=(12, 14))

    # Fig. 1a
    fig2a(
        root_folder=root_folder,
        lineage_names=lineage_names,
        dk_iar=dk_iar,
        lineage_counts=lineage_counts,
        ax=axs[0],
    )

    # Fig. 1b
    fig2b(
        proportions=proportions,
        lineage_names=lineage_names,
        dk_iar=dk_iar,
        ax=axs[1],
    )

    # Fig. 1c
    fig2c(seq_prop=seq_prop, dk_iar=dk_iar, ax=axs[2])

    # Fig. 1d
    fig2d(root_folder=root_folder, dk_iar=dk_iar, ax=axs[3])

    # Fig. 1e
    fig2e(nucleotide_diversity=nucleotide_diversity, dk_iar=dk_iar, ax=axs[4])

    # Fig. 1f
    fig2f(
        lineage_names=lineage_names,
        dk_iar=dk_iar,
        dr_list=dr_list,
        ax=axs[5],
    )

    # Formatting
    sns.despine()

    date_labs = [
        d.date().strftime("%b\n%Y")
        for d in pd.date_range(start="1/1/2021", end="1/1/2022", freq="D")
    ]

    for i in range(plt_rows):
        axs[i].label_outer()
        axs[i].xaxis.set_ticks(np.arange(N_DAYS + 1).astype(int), date_labs)

        axs[i].xaxis.set_major_locator(MonthLocator(np.arange(12).astype(int) + 1))

        if i != 1:
            axs[i].grid(alpha=0.5)
        axs[i].set_xlim([0, N_DAYS])

    # Make legend
    handles, labels = axs[1].get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    handles += [axs[0].get_legend_handles_labels()[0][-1]]
    labels += [axs[0].get_legend_handles_labels()[1][-1]]
    for i in range(2, 4):
        handles += axs[i].get_legend_handles_labels()[0]
        labels += axs[i].get_legend_handles_labels()[1]
    fig.legend(handles, labels, bbox_to_anchor=(0.75, 0.0), ncol=3, title="Legend")
    fig.tight_layout()
    plt.show()


def main():
    """
    Main function: load data and plot Figure 2
    (with extra plots on nucleotide diversity and relative growth rates)
    """
    # Folders
    args = parse_args()
    gp_folder = f"{args.root_folder}/GP_output"

    proportions, r_store, lineage_counts, dk_iar, seq_prop, nucleotide_diversity = (
        load_data(root_folder=args.root_folder, gp_folder=gp_folder)
    )

    # Get lineage names and growth rates
    lineage_names = list(r_store.columns[1:])
    lineage_filenames = [ln + "_r_boot.csv" for ln in lineage_names]
    growth_rates_boots = {}
    for i, ln in enumerate(lineage_filenames):
        growth_rates_boots[lineage_names[i]] = pd.read_csv(
            f"{gp_folder}/{ln}", index_col=0
        )

    plot_gam(seq_prop)

    plot_nucleotide_diversity(
        nucleotide_diversity, days=proportions.DaysSinceStart.values
    )

    # Part of Fig S1
    dr_list = plot_relative_growth_rates(gp_folder, lineage_names)

    fig2(
        root_folder=args.root_folder,
        proportions=proportions,
        lineage_names=lineage_names,
        lineage_counts=lineage_counts,
        dk_iar=dk_iar,
        seq_prop=seq_prop,
        nucleotide_diversity=nucleotide_diversity,
        dr_list=dr_list,
    )


if __name__ == "__main__":
    main()
