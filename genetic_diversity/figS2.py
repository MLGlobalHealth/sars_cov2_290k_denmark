"""Script to recreate Figure S2"""

# Author: Jacob Curran Sebastian (jacob.curran@sund.ku.dk)
# pylint: disable=redefined-outer, invalid-name

# Read in hamming distances, cophenetic distances and tajima d for each variant (and for all variants) - these can be calculated in Hamming.R
# Nuclueotide diversity is calculated via the function _calculate_nucleotide_diversity_single and requires a path to sequences. 

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from matplotlib.dates import MonthLocator
from joblib import delayed, Parallel
from tqdm import tqdm

from population_level_trends.data import load_data
from utils.config import N_DAYS, REF_LEN


set2 = sns.color_palette("Set2")


def parse_args():
    """Parse options to recreate Fig S2."""
    parser = ArgumentParser(
        description="Arguments for Figure S2",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-folder",
        type=str,
        help="Input folder path",
        required=True,
    )
    parser.add_argument(
        "--nd-file",
        type=str,
        help="Path to nucleotide diversity file",
        required=True,
    )
    parser.add_argument(
        "--tajima-files",
        nargs="+",
        help="list of tajima outputs from hamming.R (.csv)",
        required=True,
    )
    parser.add_argument(
        "--metadata-path",
        type=str,
        help="Metadata file",
        required=True,
    )
    parser.add_argument(
        "--coph-file",
        type=str,
        help="Path to cophenetic distance file",
        required=True,
    )
    return parser.parse_args()


def _calculate_nucleotide_diversity_single(sequences, timer=False):
    def _seqs_to_array(sequences):

        return np.array([np.array(list(str(seq))) for seq in sequences])

    sequence_array = _seqs_to_array(sequences)
    nseqs, nsites = sequence_array.shape
    nucleotides = ["A", "C", "T", "G"]
    pi_vec = np.zeros(nsites)
    ngaps = 0
    if timer == False:
        iterator = range(nsites)
    else:
        iterator = tqdm(range(nsites))
    for i in iterator:
        num_valid_sites = sum(
            [sequence_array[s, i] in nucleotides for s in range(nseqs)]
        )
        if num_valid_sites == 0:
            ngaps += 1
            continue
        pi_A = np.sum(sequence_array[:, i] == "A")
        pi_C = np.sum(sequence_array[:, i] == "C")
        pi_T = np.sum(sequence_array[:, i] == "T")
        pi_G = np.sum(sequence_array[:, i] == "G")

        summation = pi_A * (pi_C + pi_G + pi_T) + pi_C * (pi_G + pi_T) + pi_G * pi_T
        if num_valid_sites > 1:
            pi_vec[i] = summation * 2 / (num_valid_sites * (num_valid_sites - 1))

    if nsites != ngaps:
        diversity = np.sum(pi_vec) / (nsites - ngaps)

    return diversity


def calculate_nucleotide_diversity(path_to_date_file, nd_file):
    sequences_path = "./sequences.fasta"
    start = time()
    data = SeqIO.parse(sequences_path, "fasta")
    records = [(record.id, record.seq) for record in data]
    ids, sequences = list(zip(*records))
    stop = time()
    print("Data read in " + str(stop - start) + " seconds")

    sequences_by_day = []

    df = pd.read_csv(path_to_date_file, index_col=0).dropna(subset=["SampleDateTime"])

    dates = df.SampleDateTime.to_numpy()
    datetimes_first_sample = np.array(
        [datetime.strptime(date.split()[0], "%Y-%m-%d") for date in dates]
    )
    df["DateTimes"] = datetimes_first_sample
    unique_days = np.unique(datetimes_first_sample)

    for i, d in enumerate(tqdm(unique_days)):
        file_ids_by_day = df.loc[df["DateTimes"] == d].Sequence_ID.values
        sequences_by_day.append(
            [sequences[i] for i in np.where([idx in file_ids_by_day for idx in ids])[0]]
        )

    # Might be a slow step
    start = time()
    pi_by_day = Parallel(n_jobs=-2)(
        delayed(_calculate_nucleotide_diversity_single)(seqs, timer=False)
        for seqs in sequences_by_day
    )
    stop = time()
    print("Nucleotide Diversity calculated in " + str(stop - start) + " seconds")
    np.savetxt(nd_file, np.array(pi_by_day))


def get_hamming_data(metadata_path):
    metadata = pd.read_csv(metadata_path)
    hamming_dir = ""
    hamming_day = ""

    mean_hamming = np.zeros(N_DAYS)
    hamming_uq = np.zeros(N_DAYS)
    hamming_lq = np.zeros(N_DAYS)

    mean_hamming_alpha = np.zeros(N_DAYS)
    mean_hamming_delta = np.zeros(N_DAYS)
    mean_hamming_omicron = np.zeros(N_DAYS)

    for i in tqdm(range(N_DAYS)):
        hamming_df = pd.read_csv(hamming_dir + hamming_day + str(i))
        hamming_df.index = hamming_df.columns
        metadata_day = metadata[metadata["strain"].isin(hamming_df.columns.values)]
        hamming = hamming_df.values
        hamming_vals = hamming[np.tril_indices(hamming.shape[0], k=1)]
        mean_hamming[i] = np.mean(hamming_vals)
        hamming_uq[i] = np.quantile(hamming_vals, 0.95)
        hamming_lq[i] = np.quantile(hamming_vals, 0.05)

        alpha_idxs = metadata_day.loc[metadata_day["variant"] == "Alpha"].strain.values
        if len(alpha_idxs) != 0:
            hamming_alpha = (hamming_df[alpha_idxs].T[alpha_idxs]).values
            hamming_alpha_vals = hamming_alpha[
                np.tril_indices(hamming_alpha.shape[0], k=1)
            ]
            mean_hamming_alpha[i] = np.mean(hamming_alpha_vals)

        delta_idxs = metadata_day.loc[metadata_day["variant"] == "Delta"].strain.values
        if len(delta_idxs) != 0:
            hamming_delta = (hamming_df[delta_idxs].T[delta_idxs]).values
            hamming_delta_vals = hamming_delta[
                np.tril_indices(hamming_delta.shape[0], k=1)
            ]
            mean_hamming_delta[i] = np.mean(hamming_delta_vals)

        omicron_idxs = metadata_day.loc[
            metadata_day["variant"].isin(
                [
                    "Omicron (BA.1-like)",
                    "Omicron (BA.2-like)",
                    "Omicron (Unassigned)",
                    "Omicron (BA.3-like)",
                    "Omicron (BA.5-like)",
                    "Omicron (BA.4-like)",
                ]
            )
        ].strain.values
        if len(omicron_idxs) != 0:
            hamming_omicron = (hamming_df[omicron_idxs].T[omicron_idxs]).values
            hamming_omicron_vals = hamming_omicron[
                np.tril_indices(hamming_omicron.shape[0], k=1)
            ]
            mean_hamming_omicron[i] = np.mean(hamming_omicron_vals)

    return mean_hamming, mean_hamming_alpha, mean_hamming_delta, mean_hamming_omicron


def figS2a(nd_file, dk_iar, nucleotide_diversity, lw, ax):
    nd_within_lineages = pd.read_csv(nd_file, index_col=0)

    nd_alpha = nd_within_lineages["Alpha - B.1.1.7"].values
    nd_delta = nd_within_lineages["Delta - B.1.617.2"].values
    nd_omicron = nd_within_lineages["Omicron - BA.1"].values

    ax.plot(
        dk_iar.date,
        (np.array(nucleotide_diversity[1:366])),
        label="Nucleotide Diversity",
        color="tab:blue",
        linewidth=lw,
    )
    ax.set_ylabel(r"$\pi$")
    ax.plot(
        dk_iar.date[nd_alpha[:-1] > 0],
        nd_alpha[nd_alpha > 0],
        ".",
        color=set2[0],
        label="Alpha",
    )
    ax.plot(
        dk_iar.date[nd_delta[:-1] > 0],
        nd_delta[nd_delta > 0],
        ".",
        color=set2[1],
        label="Delta",
    )
    ax.plot(
        dk_iar.date[nd_omicron[:-1] > 0],
        nd_omicron[nd_omicron > 0],
        ".",
        color=set2[2],
        label="Omicron",
    )

    ax.locator_params(axis="x", nbins=12)
    # ax[3].set_ylim([0.4, 1.8])
    ax.set_title("Nucleotide Diversity")
    ax.set_xlabel("Date")


def figS2b(
    dk_iar,
    days,
    tajima_files,
    lw,
    ax,
):
    tajima_file, tajima_alpha_file, tajima_delta_file, tajima_omicron_file = (
        tajima_files
    )

    tajima_boot = pd.read_csv(tajima_file)
    tajima_boot = tajima_boot.values

    tajima_omicron = pd.read_csv(tajima_omicron_file, index_col=0).x.values
    tajima_delta = pd.read_csv(tajima_delta_file, index_col=0).x.values
    tajima_alpha = pd.read_csv(tajima_alpha_file, index_col=0).x.values

    mean_tajima = np.mean(tajima_boot[:, :], axis=0)
    uq_tajima = np.quantile(tajima_boot[:, :], 0.95, axis=0)
    lq_tajima = np.quantile(tajima_boot[:, :], 0.05, axis=0)

    ax.plot(dk_iar.date[340:], tajima_omicron[340:], ".", color=set2[2])
    ax.plot(dk_iar.date[: len(tajima_alpha)], tajima_alpha, ".", color=set2[0])
    ax.plot(dk_iar.date[91:365], tajima_delta[91:365], ".", color=set2[1])
    ax.plot(days, mean_tajima, color="tab:blue", alpha=1, linewidth=lw)
    ax.fill_between(days, lq_tajima, uq_tajima, alpha=0.1, color="tab:blue")
    ax.axhline(-2, 0, 366, linestyle="--", color="grey", alpha=0.5, linewidth=lw)
    ax.set_ylim([-4, 0])
    ax.set_ylabel("Tajima's D")
    ax.set_title("Tajima's D statistic")


def figS2c(
    metadata_path,
    dk_iar,
    lw,
    ax,
):
    mean_hamming, mean_hamming_alpha, mean_hamming_delta, mean_hamming_omicron = (
        get_hamming_data(metadata_path)
    )

    ax.plot(dk_iar.date, mean_hamming, linewidth=lw)
    alpha_plot_idxs = np.argwhere(mean_hamming_alpha > 0).flatten()
    ax.plot(
        dk_iar.date.iloc[alpha_plot_idxs],
        mean_hamming_alpha[alpha_plot_idxs],
        ".",
        color=set2[0],
    )
    delta_plot_idxs = np.argwhere(mean_hamming_delta > 0).flatten()
    ax.plot(
        dk_iar.date.iloc[delta_plot_idxs],
        mean_hamming_delta[delta_plot_idxs],
        ".",
        color=set2[1],
    )
    omicron_plot_idxs = np.argwhere(mean_hamming_omicron > 0).flatten()
    ax.plot(
        dk_iar.date.iloc[omicron_plot_idxs],
        mean_hamming_omicron[omicron_plot_idxs],
        ".",
        color=set2[2],
    )
    ax.set_ylim([0, 65])
    ax.set_title("Mean Pairwise Hamming Distance")
    ax.set_ylabel("Hamming Distance")


def figS2d(coph_file, dk_iar, lw, ax):
    cophenetic_distance = pd.read_csv(coph_file)
    mean_cophenetic = cophenetic_distance.loc[
        cophenetic_distance["region_name"] == "Full"
    ].mean_distance.values
    # copehenetic_uq = cophenetic_distance.loc[
    #     cophenetic_distance["region_name"] == "Full"
    # ].ci_high.values
    # copehenetic_lq = cophenetic_distance.loc[
    #     cophenetic_distance["region_name"] == "Full"
    # ].ci_low.values

    ax.plot(dk_iar.date, mean_cophenetic * REF_LEN, color="tab:blue", linewidth=lw)
    ax.set_title("Mean Pairwise Cophenetic Distance")
    ax.set_ylabel("Cophenetic Distance")
    ax.set_xlabel("Date")
    ax.set_ylim([0, 65])


def figS2(
    dk_iar,
    nucleotide_diversity,
    days,
    nd_file,
    tajima_files,
    metadata_path,
    coph_file,
    start_date="1/1/2021",
    end_date="1/1/2022",
):
    date_labs = [
        d.date().strftime("%b\n%Y")
        for d in pd.date_range(start=start_date, end=end_date, freq="D")
    ]

    diversity_output_file = ""

    plt.style.use("nature")
    font_size = 10
    plt.rcParams["font.size"] = font_size
    plt.rcParams["legend.title_fontsize"] = font_size
    plt.rcParams["axes.labelsize"] = font_size
    plt.rcParams["legend.fontsize"] = font_size
    plt.rcParams["xtick.labelsize"] = font_size
    plt.rcParams["xtick.labelsize"] = font_size

    lw = 1.5

    plt_rows = 4
    fig, axs = plt.subplots(nrows=plt_rows, ncols=1, figsize=(10, 15))

    figS2a(
        nd_file=nd_file,
        dk_iar=dk_iar,
        nucleotide_diversity=nucleotide_diversity,
        lw=lw,
        ax=axs[0],
    )

    figS2b(dk_iar=dk_iar, days=days, tajima_files=tajima_files, lw=lw, ax=axs[1])

    figS2c(metadata_path=metadata_path, dk_iar=dk_iar, lw=lw, ax=axs[2])

    figS2d(coph_file=coph_file, dk_iar=dk_iar, lw=lw, ax=axs[3])

    handles, labels = axs[0].get_legend_handles_labels()
    handles = handles[1:]
    labels = labels[1:]
    lgd = fig.legend(handles, labels, bbox_to_anchor=(1.1, 0.97), title="Variant")

    for i in range(plt_rows):
        axs[i].label_outer()
        axs[i].xaxis.set_ticks(np.arange(366).astype(int), date_labs)

        axs[i].xaxis.set_major_locator(MonthLocator(np.arange(12).astype(int) + 1))

        axs[i].grid(alpha=0.5)
        axs[i].set_xlim([0, 365])

    plt.tight_layout()
    plt.savefig(diversity_output_file, bbox_extra_artists=(lgd,), bbox_inches="tight")


def main():
    """
    Main function: load data and plot Figure S2
    """
    # Folders
    args = parse_args()
    gp_folder = f"{args.root_folder}/GP_output"

    proportions, _, _, dk_iar, _, nucleotide_diversity = load_data(
        root_folder=args.root_folder, gp_folder=gp_folder
    )

    days = proportions.DaysSinceStart.values
    figS2(
        dk_iar=dk_iar,
        nucleotide_diversity=nucleotide_diversity,
        days=days,
        nd_file=args.nd_file,
        tajima_files=args.tajima_files,
        metadata_path=args.metadata_path,
        coph_file=args.coph_file,
    )


if __name__ == "__main__":
    main()
