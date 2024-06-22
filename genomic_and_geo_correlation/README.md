# Figure 6: relationship between geographic and genomic distances alongside mean cophenetic distances within and between households, and mean cophenetic distance over time by region

## OSRM distance matrix extraction

### Dependencies (Julia)

* ArgParse
* CSV
* Combinatorics
* DataFrames
* HTTP
* JLD2
* JSON
* ProgressMeter

### Option 1: use Apptainer to run the OSRM API

* Use Apptainer to run the OSRM API (e.g., on port 5000)
* Once port 5000 is open and listening, we can run build_distance_matrix.jl

```bash
export JULIA_NUM_THREADS=16 # change number as desired

julia build_distance_matrix.jl
```

### Option 2: use downloaded maps

#### Creating maps

* Install ```osmosis```: <https://wiki.openstreetmap.org/wiki/Osmosis/Installation#Linux>

* Download the latest national map from: <https://download.geofabrik.de/europe/denmark-latest.osm.pbf> and save it under ```data/maps```

```bash
osmosis --read-pbf data/maps/denmark-latest.osm.pbf --bb left=8.189517 bottom=56.534547 right=11.225991 top=57.760255 --write-pdf data/maps/nordjylland.osm.pbf
osmosis --read-pbf data/maps/denmark-latest.osm.pbf --bb left=8.078876 bottom=55.644379 right=11.664191 top=56.843257 --write-pbf data/maps/midtjylland.osm.pbf
osmosis --read-pbf data/maps/denmark-latest.osm.pbf --bb left=8.063203 bottom=54.718281 right=10.995552 top=55.953250 --write-pbf data/maps/syddanmark.osm.pbf
osmosis --read-pbf data/maps/denmark-latest.osm.pbf --bb left=10.814805 bottom=54.544406 right=12.645516 top=56.017306 --write-pbf data/maps/sjælland.osm.pbf
```

#### Script

```bash
export JULIA_NUM_THREADS=16 # change number as desired

julia build_distance_matrix_v2.jl
```

## Downstream analyses (R)

### Dependencies (R)

```R
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, ggplotify,
  phangorn, heatmaply, RColorBrewer, graphics, purrr, future.apply, geosphere, patchwork, coefplot,
  adephylo, biglm, pheatmap, hdf5r, JuliaCall
)
```
