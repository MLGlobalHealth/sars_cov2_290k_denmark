# Data for Molecular Change Analysis

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Reading relevant packages ----------
library(ape)

# Reading in data
sequenced_individual_detailed_metadata <- readRDS("data/sequenced_individual_detailed_metadata.RDS")
distance_tree <- read.tree("data/tree/distance_tree.tree")
time_tree <- read.tree("data/tree/time_tree.tree")

# Checking relationship between age and branch lengths -------
# Extract tip labels from distance tree
tip_labels <- distance_tree$tip.label
corresponding_PERSON_ID <- sequenced_individual_detailed_metadata$PERSON_ID[
  match(tip_labels, sequenced_individual_detailed_metadata$strain)
]

filtered_tree <- distance_tree
filtered_tree$tip.label <- corresponding_PERSON_ID

ages <- sequenced_individual_detailed_metadata$age_at_infection[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

regions <- sequenced_individual_detailed_metadata$REGIONSKODE[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

sex <- sequenced_individual_detailed_metadata$KOEN[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

partial_vacc <- sequenced_individual_detailed_metadata$Partially_vaccinated_at_exposure[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

fully_vacc <- sequenced_individual_detailed_metadata$Fully_vaccinated_at_exposure[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]
date <- sequenced_individual_detailed_metadata$date[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

major_variant <- sequenced_individual_detailed_metadata$variant[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

PERSON_ID <- sequenced_individual_detailed_metadata$PERSON_ID[
  match(corresponding_PERSON_ID, sequenced_individual_detailed_metadata$PERSON_ID)
]

# Keep only the branch lengths going to the tips and call them edge_lengths
edge_lengths <- setNames(
  filtered_tree$edge.length[sapply(seq_along(filtered_tree$tip.label),
    function(x, y) which(y == x),
    y = filtered_tree$edge[, 2]
  )],
  filtered_tree$tip.label
)

# Creating one dataframe with all lengths, one with just non-zero branch lengths -------
# All Data
analysis_data <- data.frame(
  PERSON_ID = PERSON_ID, ages = ages, regions = regions,
  sex = sex, date = date, partial_vacc = partial_vacc,
  fully_vacc = fully_vacc,
  branch_lengths = edge_lengths,
  major_variant = major_variant
)
analysis_data$month <- month(analysis_data$date)


# Adding tip lengths from time tree ------
# Extract tip labels from time tree
time_tip_labels <- time_tree$tip.label
corresponding_PERSON_ID_time <- sequenced_individual_detailed_metadata$PERSON_ID[
  match(time_tip_labels, sequenced_individual_detailed_metadata$strain)
]
filtered_time_tree <- time_tree
filtered_time_tree$tip.label <- corresponding_PERSON_ID_time

# Keep only the branch lengths going to the tips and call them edge_lengths_time
edge_lengths_time <- setNames(
  filtered_time_tree$edge.length[sapply(seq_along(filtered_time_tree$tip.label),
    function(x, y) which(y == x),
    y = filtered_time_tree$edge[, 2]
  )],
  filtered_time_tree$tip.label
)

# Make sure the tips correspond to the right distance_tree tips
edge_lengths_time <- edge_lengths_time[match(analysis_data$PERSON_ID, names(edge_lengths_time))]

# Now, add the branch lengths to the analysis_data dataframe
analysis_data$branch_lengths_time <- edge_lengths_time
analysis_data$rate <- (analysis_data$branch_lengths) / (analysis_data$branch_lengths_time)
write.csv(analysis_data, file = "data/sequenced_individual_detailed_metadata.csv")
