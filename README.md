# High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark

Code for "High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark"

Link to paper: <https://doi.org/10.1038/s41467-024-51371-0>

## Overview

The presented code aims to reproduce the analyses from "High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark". In brief, the goal of the code is to identify population-level trends (e.g., case counts, nucleotide diversity), create phylogenetic trees (e.g., for all sequences and for each variant-specific clade) and analyze the relationship between various demographic characteristics and molecular change (e.g., examining differences in tip lengths and rates between different demographic groups) using 290k sequences from Denmark in 2021.

## Folders

* Figure 1: workflow (see manuscript)
* Figure 2, S12: ```population_level_trends```
* Figure 3: ```phylogenetic tree```
* Figures 4, S4, S5, S6, S7: ```clade_characterization```
* Figure 5, S8, S9, S10: ```evolutionary_rates```
* Figure 6, S11: ```genomic_and_geo_correlation```
* Figure S1: ```growth_rates```
* Figure S2, S3: ```genetic_diversity```

To reproduce everything, install the .yml with conda/mamba:

```bash
conda env create -f environment.yml
```

To reproduce certain figures individually, consult the README files in each corresponding folder.

Full newick files for each clade and for the whole tree are found in ```phylogenetic_newick_trees```

## General Data

Some data is available under the ```data``` folder. Due to confidentiality restrictions, dummy/synthetic data has been generated; they do not correspond to any real data, but were generated to allow for testing of the code's functionality.

* ```data/BEAST_XML_files``` contains the XML files necessary to run BEAST; for confidentiality reasons, the XML files do not include taxon or sequence data
* ```clade_characterization/data/phylogenetic_newick_trees``` contains phylogenetic trees (in Newick format) where the sequence IDs are anonymized
* ```data/other``` contains other publicly available data

## Script-Specific Data

Synthetic data is available under the ```data``` folder in each of the main sub-folders. Due to confidentiality restrictions, dummy/synthetic data to test the code has been generated instead.


### Genomes

259,106 high-quality SARS-CoV-2 consensus genomes used in this study are available on the GISAID’s EpiCoV database under the EPI-SET accession number EPI_SET_240423qn.

## Software overview

R (4.2.3), Python (3.10.10), Julia, BEAST (1.10.4) for clade-specific phylogenetic tree inference, MAPLE (0.3.1), Chronumental (0.0.60) for creating time trees from the MAPLE distance tree, MAFFT (7.520) for sequence alignment

## Authors

* Mark Khurana (<mark.khurana@sund.ku.dk>)
* Jacob-Curran Sebastian (<jacob.curran@sund.ku.dk>)
* Neil Scheidwasser (<neil.clow@sund.ku.dk>)

## License

Apache 2.0 License

## Citation

Please cite the paper as:

```tex
@article{khurana2024high,
  title={High-resolution epidemiological landscape from\~{} 290,000 SARS-CoV-2 genomes from Denmark},
  author={Khurana, Mark P and Curran-Sebastian, Jacob and Scheidwasser, Neil and Morgenstern, Christian and Rasmussen, Morten and Fonager, Jannik and Stegger, Marc and Tang, Man-Hung Eric and Juul, Jonas L and Escobar-Herrera, Leandro Andr{\'e}s and others},
  journal={Nature Communications},
  volume={15},
  number={1},
  pages={7123},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```
