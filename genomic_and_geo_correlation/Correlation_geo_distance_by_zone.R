# Correlation, geographic region with cophenetic distance by zone (urban or rural)

# Reading packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr,
  lubridate, rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica,
  DSTdb, DBI, parallel, ggsignif, Rcpp, geosphere, biglm, graphics,
  pheatmap, viridis, patchwork, coefplot, ggpubr
)


# Loading basic data -----
sequenced_individuals <- readRDS(file = "")
final_tree <- read.tree("")
final_distance_tree <- read.tree("")

# Filter sequenced_individuals dataframe by zone (urban or land/countryside) -----
byzone_data <- sequenced_individuals[sequenced_individuals$ZONE == "Byzone", ]
landzone_data <- sequenced_individuals[sequenced_individuals$ZONE == "Landzone", ]


# Functions for later --------
# Create a function to calculate Haversine distance for a pair of coordinates
calculate_haversine_distance <- function(coord1, coord2) {
  distHaversine(coord1, coord2, r = 6378137)
}

# Create a function to calculate the absolute time difference in days between two dates
calculate_time_difference <- function(date1, date2) {
  as.numeric(difftime(date1, date2, units = "days"))
}


# BYZONE, 18111 people -------
# Extract PERSON_ID values for the sampled rows
sampled_person_ids <- as.character(byzone_data$strain)
matched_strains <- match(final_distance_tree$tip.label, sampled_person_ids)
tips_to_keep_byzone <- final_distance_tree$tip.label[!is.na(matched_strains)]
tips_to_keep_byzone <- as.character(tips_to_keep_byzone)
pruned_byzone_tree <- keep.tip(final_distance_tree, tips_to_keep_byzone)
cophenetic_distances_byzone <- cophenetic(pruned_byzone_tree)

# Get the order of names from the cophenetic_distances matrix
order_of_names_byzone <- rownames(cophenetic_distances_byzone)
sampled_data_reordered <- byzone_data[match(order_of_names_byzone, byzone_data$strain), ]

# Distance matrix for geographic locations
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
coordinates <- sampled_data_reordered[, c("longitude", "latitude")]
geographic_distance_matrix_byzone <- foreach(
  i = seq_len(nrow(sampled_data_reordered)),
  .combine = rbind
) %dopar% {
  # Load the geosphere library inside the parallel block
  library(geosphere)

  sapply(seq_len(nrow(sampled_data_reordered)), function(j) {
    calculate_haversine_distance(coordinates[i, ], coordinates[j, ])
  })
}
stopCluster(cl)

# Time between people
cl <- makeCluster(num_cores)
registerDoParallel(cl)
time_distance_matrix_byzone <- foreach(i = seq_len(nrow(sampled_data_reordered)), .combine = rbind) %dopar% {
  sapply(seq_len(nrow(sampled_data_reordered)), function(j) {
    calculate_time_difference(sampled_data_reordered$date[i], sampled_data_reordered$date[j])
  })
}
stopCluster(cl)

saveRDS(sampled_data_reordered, file = "")
saveRDS(cophenetic_distances_byzone, file = "")
saveRDS(geographic_distance_matrix_byzone, file = "")
saveRDS(time_distance_matrix_byzone, file = "")






# LANDZONE, 1817 people -------
# Extract PERSON_ID values for the sampled rows
sampled_person_ids <- as.character(landzone_data$strain)
matched_strains <- match(final_distance_tree$tip.label, sampled_person_ids)
tips_to_keep_landzone <- final_distance_tree$tip.label[!is.na(matched_strains)]
tips_to_keep_landzone <- as.character(tips_to_keep_landzone)
pruned_landzone_tree <- keep.tip(final_distance_tree, tips_to_keep_landzone)
cophenetic_distances_landzone <- cophenetic(pruned_landzone_tree)

# Get the order of names from the cophenetic_distances matrix
order_of_names_landzone <- rownames(cophenetic_distances_landzone)
sampled_data_reordered_landzone <- landzone_data[match(order_of_names_landzone, landzone_data$strain), ]

# Distance matrix for geographic locations
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
coordinates <- sampled_data_reordered_landzone[, c("longitude", "latitude")]
geographic_distance_matrix_landzone <- foreach(i = seq_len(nrow(sampled_data_reordered_landzone)), .combine = rbind) %dopar% {
  library(geosphere)

  sapply(seq_len(nrow(sampled_data_reordered_landzone)), function(j) {
    calculate_haversine_distance(coordinates[i, ], coordinates[j, ])
  })
}
stopCluster(cl)

# Time between people
cl <- makeCluster(num_cores)
registerDoParallel(cl)
time_distance_matrix_landzone <- foreach(i = seq_len(nrow(sampled_data_reordered_landzone)), .combine = rbind) %dopar% {
  sapply(seq_len(row(sampled_data_reordered_landzone)), function(j) {
    calculate_time_difference(
      sampled_data_reordered_landzone$date[i],
      sampled_data_reordered_landzone$date[j]
    )
  })
}
stopCluster(cl)

saveRDS(sampled_data_reordered_landzone, file = "")
saveRDS(cophenetic_distances_landzone, file = "")
saveRDS(geographic_distance_matrix_landzone, file = "")
saveRDS(time_distance_matrix_landzone, file = "")





# Postcode data for major cities ------
postcodes <- read.csv(file = "")
copenhagen_area <- subset(postcodes, grepl("^København|^Frederiksberg", Commune))
copenhagen_postcodes <- unique(copenhagen_area$Postcode)

aarhus_area <- subset(postcodes, grepl("^Aarhus", Commune))
aarhus_postcodes <- unique(aarhus_area$Postcode)

odense_area <- subset(postcodes, grepl("^Odense", Commune))
odense_postcodes <- unique(odense_area$Postcode)

copenhagen <- data.frame(
  city = "Copenhagen",
  postcode = paste(copenhagen_postcodes, collapse = ", ")
)
aarhus <- data.frame(city = "Aarhus", postcode = paste(aarhus_postcodes, collapse = ", "))
odense <- data.frame(city = "Odense", postcode = paste(odense_postcodes, collapse = ", "))
# Combining data frames
cities <- rbind(copenhagen, aarhus, odense)
write.csv(cities, file = "", row.names = FALSE)




# Only Copenhagen, 3416 individuals --------
copenhagen_data <- sequenced_individuals[sequenced_individuals$POSTNR %in% copenhagen_postcodes, ]
# Extract PERSON_ID values for the sampled rows
sampled_person_ids <- as.character(copenhagen_data$strain)
matched_strains <- match(final_distance_tree$tip.label, sampled_person_ids)
tips_to_keep_copenhagen <- final_distance_tree$tip.label[!is.na(matched_strains)]
tips_to_keep_copenhagen <- as.character(tips_to_keep_copenhagen)
pruned_copenhagen_tree <- keep.tip(final_distance_tree, tips_to_keep_copenhagen)
cophenetic_distances_copenhagen <- cophenetic(pruned_copenhagen_tree)

# Get the order of names from the cophenetic_distances matrix
order_of_names_copenhagen <- rownames(cophenetic_distances_copenhagen)
sampled_data_reordered_copenhagen <- copenhagen_data[match(order_of_names_copenhagen, copenhagen_data$strain), ]

# Distance matrix for geographic locations
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
coordinates <- sampled_data_reordered_copenhagen[, c("longitude", "latitude")]
geographic_distance_matrix_copenhagen <- foreach(
  i = seq_len(nrow(sampled_data_reordered_copenhagen)),
  .combine = rbind
) %dopar% {
  library(geosphere)

  sapply(seq_len(nrow(sampled_data_reordered_copenhagen)), function(j) {
    calculate_haversine_distance(coordinates[i, ], coordinates[j, ])
  })
}
stopCluster(cl)

# Time between people
cl <- makeCluster(num_cores)
registerDoParallel(cl)
time_distance_matrix_copenhagen <- foreach(
  i = seq_len(nrow(sampled_data_reordered_copenhagen)),
  .combine = rbind
) %dopar% {
  sapply(seq_len(nrow(sampled_data_reordered_copenhagen)), function(j) {
    calculate_time_difference(
      sampled_data_reordered_copenhagen$date[i],
      sampled_data_reordered_copenhagen$date[j]
    )
  })
}
stopCluster(cl)

saveRDS(sampled_data_reordered_copenhagen, file = "")
saveRDS(cophenetic_distances_copenhagen, file = "")
saveRDS(geographic_distance_matrix_copenhagen, file = "")
saveRDS(time_distance_matrix_copenhagen, file = "")
