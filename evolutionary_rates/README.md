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

One can also use the synthpop package to generate synthetic data and run the ```synthesise.R``` script:

```r
install.packages(synthpop)

source("synthesise.R")
```
