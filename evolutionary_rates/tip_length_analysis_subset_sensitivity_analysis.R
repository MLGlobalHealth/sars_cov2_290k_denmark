# Data for Molecular Change Analysis, subset tree from April to November

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Load packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table,
  dbplyr, lubridate, rlang, foreach, doParallel, parallel, ggsignif,
  caper, picante, mgcv, patchwork, coefplot, ggpubr, stargazer, lme4
)


# Reading in data
sequenced_individual_detailed_metadata <- readRDS(file = "data/sequenced_individual_detailed_metadata.RDS")
distance_tree <- read.tree("data/tree/distance_tree.tree")
time_tree <- read.tree("data/tree/time_tree.tree")

# Filter metadata to include only sequences between April and November
filtered_metadata <- sequenced_individual_detailed_metadata %>%
  filter(date >= as.Date("2021-04-01") & date < as.Date("2021-11-01"))

# Get the IDs of the filtered individuals
filtered_ids <- filtered_metadata$strain

# Filter the tips from the distance tree to include only those individuals
filtered_distance_tree <- keep.tip(distance_tree, tip = filtered_ids)

# Filter the tips from the time tree to include only those individuals
filtered_time_tree <- keep.tip(time_tree, tip = filtered_ids)

# Checking relationship between age and branch lengths -------
# Extract tip labels from distance tree
tip_labels <- filtered_distance_tree$tip.label
corresponding_PERSON_ID <- filtered_metadata$PERSON_ID[match(tip_labels, filtered_metadata$strain)]
filtered_tree <- filtered_distance_tree
filtered_tree$tip.label <- corresponding_PERSON_ID

# Retrieve other metadata corresponding to the filtered tip labels
ages <- filtered_metadata$age_at_infection[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
regions <- filtered_metadata$REGIONSKODE[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
sex <- filtered_metadata$KOEN[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
partial_vacc <- filtered_metadata$Partially_vaccinated_at_exposure[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
fully_vacc <- filtered_metadata$Fully_vaccinated_at_exposure[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
date <- filtered_metadata$date[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
major_variant <- filtered_metadata$variant[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]
PERSON_ID <- filtered_metadata$PERSON_ID[
  match(corresponding_PERSON_ID, filtered_metadata$PERSON_ID)
]


# Keep only the branch lengths going to the tips and call them edge_lengths
edge_lengths <- setNames(filtered_tree$edge.length[sapply(1:length(filtered_tree$tip.label), function(x, y) which(y == x), y = filtered_tree$edge[, 2])], filtered_tree$tip.label)

# Create dataframe for tips and features
analysis_data <- data.frame(
  PERSON_ID = PERSON_ID, ages = ages, regions = regions,
  sex = sex, date = date, partial_vacc = partial_vacc,
  fully_vacc = fully_vacc,
  branch_lengths = edge_lengths,
  major_variant = major_variant
)
analysis_data$month <- month(analysis_data$date)

# Extract tip labels from time tree
time_tip_labels <- filtered_time_tree$tip.label
corresponding_PERSON_ID_time <- filtered_metadata$PERSON_ID[match(time_tip_labels, filtered_metadata$strain)]
filtered_time_tree$tip.label <- corresponding_PERSON_ID_time
# Keep only the branch lengths going to the tips and call them edge_lengths_time
edge_lengths_time <- setNames(filtered_time_tree$edge.length[sapply(1:length(filtered_time_tree$tip.label), function(x, y) which(y == x), y = filtered_time_tree$edge[, 2])], filtered_time_tree$tip.label)
# Make sure the tips correspond to the right distance_tree tips
edge_lengths_time <- edge_lengths_time[match(analysis_data$PERSON_ID, names(edge_lengths_time))]
# Now, add the branch lengths to the analysis_data dataframe
analysis_data$branch_lengths_time <- edge_lengths_time
analysis_data$rate <- (analysis_data$branch_lengths) / (analysis_data$branch_lengths_time)

write.csv(analysis_data, file = "")
