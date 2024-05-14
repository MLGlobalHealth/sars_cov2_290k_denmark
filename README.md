# SARS_CoV2_290k_Denmark

Code for "High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark" (under review; available soon)

## Overview
The code here aims to reproduce the analyses from "High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark". In brief, the goal of the code is to identify population-level trends (e.g., case counts, nucleotide diversity), create phylogenetic trees (e.g., for all sequences and for each variant-specific clade) and analyze the relationship between various demographic characteristics and molecular change (e.g., examining differences in tip lengths and rates between different demographic groups) using 290k sequences from Denmark in 2021.

## Folders

* Figure 1: workflow (see manuscript)
* Figure 2: ```population_level_trends```
* Figure 3: ```phylogenetic tree```
* Figures 4, 5, S4, S5, S7: ```clade_characterization```
* Figure 6, S8, S9: ```evolutionary_rates```
* Figure 7, S10: ```genomic_and_geo_correlation```
* Figure S1: ```growth_rates```
* Figure S2, S3: ```genetic_diversity```

To reproduce everything, install the .yml with conda/mamba:

```bash
conda env create -f environment.yml
```

To reproduce certain figures individually, consult the README files in each corresponding folder.

Full newick files for each clade and for the whole tree are found in ```phylogenetic_newick_trees```

## Data availability
259,106 high-quality SARS-CoV-2 consensus genomes used in this study are available on the GISAID’s EpiCoV database under the EPI-SET accession number EPI_SET_240423qn.

## Software overview
R (4.2.3)
Python (3.10.10)
Julia
BEAST (1.10.4) for clade-specific phylogenetic tree inference
MAPLE (0.3.1)
Chronumental (0.0.60) for creating time trees from the MAPLE distance tree
MAFFT (7.520) for sequence alignment

## Authors

* Mark Khurana (<mark.khurana@sund.ku.dk>)
* Jacob-Curran Sebastian (<jacob.curran@sund.ku.dk>)
* Neil Scheidwasser (<neil.clow@sund.ku.dk>)

## License
Apache 2.0 License
