# Data for correlation for household subset

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, ggplotify,
  phangorn, heatmaply, RColorBrewer, graphics, purrr, future.apply, geosphere, patchwork, coefplot,
  adephylo, biglm, pheatmap
)


# Loading basic data -----
sequenced_individuals <- readRDS(file = "data/sequenced_individuals.RDS")
all_individual <- readRDS(file = "data/all_individual.RDS")
final_tree <- read.tree("data/tree/final_tree.tree")
final_distance_tree <- read.tree("data/tree/final_distance_tree.tree")


# Parallelized method, 3 people in a household ----------

# Filter rows where there are THREE people in the household, for computational feasiblity
filtered_household_infection_sets_three_people <- sequenced_individuals %>%
  group_by(UNIQUE_ADDRESS_ID) %>%
  filter(n_distinct(PERSON_ID) == 3 & n() == 3) %>%
  ungroup()
# Select 1000 unique household IDs
random_households <- sample(unique(filtered_household_infection_sets_three_people$UNIQUE_ADDRESS_ID), 1000, replace = FALSE)
filtered_household_infection_sets_three_people <- filtered_household_infection_sets_three_people %>%
  filter(UNIQUE_ADDRESS_ID %in% random_households) %>%
  ungroup()
saveRDS(filtered_household_infection_sets_three_people, file = "filtered_household_infection_sets_three_people.RDS")

# Now we have 1000 households with exactly 3 people testing positive in 2021

# Find all unique people in all households
all_people <- unique(filtered_household_infection_sets_three_people$strain)

# Use keep.tip to subset the tree to include only these people
subtree <- keep.tip(final_distance_tree, tip = all_people)

# Get unique households with 2-3 people
unique_households <- unique(filtered_household_infection_sets_three_people$UNIQUE_ADDRESS_ID)

# Calculate time differences:
dates <- filtered_household_infection_sets_three_people$date
date_differences_matrix <- matrix(NA, nrow = length(all_people), ncol = length(all_people), dimnames = list(all_people, all_people))
for (i in 1:length(all_people)) {
  for (j in 1:length(all_people)) {
    if (j >= i) {
      person1 <- all_people[i]
      person2 <- all_people[j]

      # Extract dates for the two individuals
      date1 <- dates[all_people == person1]
      date2 <- dates[all_people == person2]

      # Calculate the date difference and store it in the matrix
      date_difference <- as.numeric(abs(difftime(date2, date1, units = "days")))
      date_differences_matrix[person1, person2] <- date_difference
      date_differences_matrix[person2, person1] <- date_difference # Negative for person2 to person1
    }
  }
}

# Calculate cophenetic distances
cophenetic_distances <- cophenetic(subtree)

# Preallocate distance matrices
n <- length(unique_households)
distance_matrix_average_three_people <- matrix(NA, n, n, dimnames = list(unique_households, unique_households))
distance_matrix_time_three_people <- matrix(NA, n, n, dimnames = list(unique_households, unique_households))

# Iterate over unique households
for (i in 1:n) {
  for (j in 1:n) {
    if (j <= i) {
      # Extract relevant strains for people in the ith and jth households
      strains_i <- filtered_household_infection_sets_three_people$strain[filtered_household_infection_sets_three_people$UNIQUE_ADDRESS_ID == unique_households[i]]
      strains_j <- filtered_household_infection_sets_three_people$strain[filtered_household_infection_sets_three_people$UNIQUE_ADDRESS_ID == unique_households[j]]

      indices_i <- which(rownames(cophenetic_distances) %in% strains_i)
      indices_j <- which(rownames(cophenetic_distances) %in% strains_j)

      # Subset the cophenetic matrix to include only rows and columns with combined strains
      distances <- cophenetic_distances[indices_i, indices_j]
      times <- date_differences_matrix[indices_i, indices_j]

      mean_distance <- mean(distances, na.rm = TRUE)
      time_average <- mean(times, na.rm = TRUE)

      # Place mean and sum in the distance matrices
      distance_matrix_average_three_people[i, j] <- mean_distance
      distance_matrix_time_three_people[i, j] <- time_average
    }
  }
}

# Save distance matrices
# saveRDS(distance_matrix_average_three_people, file = "")
# saveRDS(distance_matrix_time_three_people, file = "")
