# Figure 5: Analysis of substitution rate variability

## Dependencies

### Installation (Python)

```python
pip install -r requirements.txt
```

### Installation (R)

For analysis:

```r
install.packages(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table,
  dbplyr, lubridate, rlang, foreach, doParallel, parallel,
  ggsignif, caper, picante, mgcv, patchwork, coefplot, ggpubr,
  stargazer, lme4
)
```

### Data

Due to privacy and confidentiality restrictions, we currently cannot provide the true dataset, even de-identified. Thus, we provide "dummy" data to show the code can be executed.

One can also use the synthpop package to generate synthetic data and run the ```synthesise.R``` script:

```r
install.packages(synthpop)

source("synthesise.R")
```
