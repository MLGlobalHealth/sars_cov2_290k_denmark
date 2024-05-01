"""Data loading."""

import numpy as np
import pandas as pd


def load_data(root_folder, gp_folder):
    """Load the data necessary for Figure 2

    Returns
    -------
    proportions : pandas.DataFrame
        Variant proportions
    r_store : pandas.DataFrame
        Growth rates
    lineage_counts : pandas.DataFrame
        Daily counts of each lineage
    dk_iar : pandas.DataFrame
        Denmark infection ascertainment rates (IARs)
    seq_prop : numpy.ndarray
        Proportion of PCR positives
    nucleotide_diversity : numpy.ndarray
        Nucleotide diversity
    """
    # Variant proportions
    probabilities_file = f"{gp_folder}/Probabilities_new.csv"
    proportions = pd.read_csv(probabilities_file)
    # Drop first row and reorder columns
    proportions.drop(0, inplace=True)
    proportions = proportions.iloc[:, [0, 9, 1, 2, 3, 4, 5, 6, 7, 8]]

    # Growth rates (and lineage names)
    r_store = pd.read_csv(f"{gp_folder}/r_store.csv")

    # Lineage counts
    lineage_counts = pd.read_csv(f"{gp_folder}/day_counts_new.csv").drop(0)
    lineage_counts[["All", "Alpha - B.1.1.7", "Other"]].head()

    # DK IARs
    dk_iar = pd.read_csv(f"{root_folder}/denmark_rt_iar.csv")
    dk_iar.query("date >= '2021-01-01'", inplace=True)

    # Proportion of PCR positives
    seq_prop = np.genfromtxt(f"{root_folder}/proportion_of_pcrpositives_sequenced.csv")

    # Nucleotide diversity
    nucleotide_diversity = np.genfromtxt(
        f"{root_folder}/GP_output/nucleotide_diversity.csv"
    )

    return proportions, r_store, lineage_counts, dk_iar, seq_prop, nucleotide_diversity
