# sars_cov2_290k_denmark

Code for "High-Resolution Epidemiological Landscape from 290K SARS-CoV-2 Genomes from Denmark" (under review; available soon)

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

## Authors

* Mark Khurana (<mark.khurana@sund.ku.dk>)
* Jacob-Curran Sebastian (<jacob.curran@sund.ku.dk>)
* Neil Scheidwasser (<neil.clow@sund.ku.dk>)
