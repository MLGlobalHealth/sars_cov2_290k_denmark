# Data generation for correlation between genetic and geographic distance

# Importing packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, ggplotify,
  phangorn, heatmaply, RColorBrewer, graphics, purrr, future.apply, geosphere, patchwork, coefplot,
  adephylo, biglm, pheatmap, osrm, oce, rgdal, sp, osmextract
)


# Data Preparation (Tree and Metadata) ------------------------------
sequenced_individuals <- readRDS(file = "")
all_individual <- read.csv(file = "")
distance_tree <- read.tree("")

# Change tip labels to PERSON_ID values
strain_person_id_mapping <- setNames(all_individual$PERSON_ID, all_individual$strain)
matched_person_ids <- strain_person_id_mapping[distance_tree$tip.label]
distance_tree$tip.label <- matched_person_ids

# Cophenetic distance matrix for 20,000 random samples ----------
# Take 20,000 random rows from sequenced_individuals
set.seed(123)
# Initialize variables
sampled_data <- data.frame()
sampled_person_ids <- character(0)
# Loop until you have sampled 20,000 unique rows
while (nrow(sampled_data) < 20000) {
  # Sample a new row
  new_row <- sequenced_individuals[sample(nrow(sequenced_individuals), 1), ]

  # Check if the PERSON_ID is already in the sampled data
  if (!(new_row$PERSON_ID %in% sampled_person_ids)) {
    # If not, add the row to the sampled data and PERSON_ID to the vector
    sampled_data <- rbind(sampled_data, new_row)
    sampled_person_ids <- c(sampled_person_ids, new_row$PERSON_ID)
  }
}

# Extract PERSON_ID values for the sampled rows
sampled_person_ids <- sampled_data$PERSON_ID
sampled_person_ids <- as.character(sampled_data$PERSON_ID)
# Match the 'PERSON_ID' values in the tree with the corresponding 'strain' values
matched_strains <- match(distance_tree$tip.label, sampled_person_ids)
# Identify tips to keep based on matched_strains
tips_to_keep <- distance_tree$tip.label[!is.na(matched_strains)]
tips_to_keep <- as.character(tips_to_keep)
# Prune the tree to retain only the sampled individuals
pruned_tree <- keep.tip(distance_tree, tips_to_keep)
# Calculate the cophenetic distance matrix
cophenetic_distances <- cophenetic(pruned_tree)

# Re-ordering the sampled data to match the cophenetic distance matrix ------
# Get the order of names from the cophenetic_distances matrix
order_of_names <- rownames(cophenetic_distances)
sampled_data_reordered <- sampled_data[match(order_of_names, sampled_data$PERSON_ID), ]
saveRDS(sampled_data_reordered, file = "sampled_subset_data.rds")
write.csv(sampled_data_reordered, file = "sampled_subset_data.csv")
sampled_data <- sampled_data_reordered

# Converting UTM/ETRS to longitude and latitude -------
sampled_data$longitude_and_latitude <- paste(sampled_data$ETRS89_EAST.x, sampled_data$ETRS89_NORTH.x, sep = ", ")
# Convert UTM coordinates to longitude and latitude
coordinates <- strsplit(as.character(sampled_data$longitude_and_latitude), ", ")
coordinates <- matrix(unlist(coordinates), ncol = 2, byrow = TRUE)
# Extract east and north values
easting <- as.numeric(coordinates[, 1])
northing <- as.numeric(coordinates[, 2])
# Convert UTM to lonlat using oce's utm2lonlat function
lonlat <- utm2lonlat(easting = easting, northing = northing, zone = 32, hemisphere = "N")
# Add the lonlat values as new columns
sampled_data$longitude <- lonlat$longitude
sampled_data$latitude <- lonlat$latitude
# month
sampled_data$date <- as.Date(sampled_data$date) # Ensure 'date' is in Date format
sampled_data$month <- month(sampled_data$date)
head(sampled_data)
# saveRDS(sampled_data, file="sampled_subset_data.rds")
# write.csv(sampled_data, file="sampled_subset_data.csv")


# Converting UTM/ETRS to longitude and latitude for the big metadata set  -------
sequenced_individuals$longitude_and_latitude <- paste(sequenced_individuals$ETRS89_EAST.x, sequenced_individuals$ETRS89_NORTH.x, sep = ", ")
# Convert UTM coordinates to longitude and latitude
coordinates <- strsplit(as.character(sequenced_individuals$longitude_and_latitude), ", ")
coordinates <- matrix(unlist(coordinates), ncol = 2, byrow = TRUE)
# Extract easting and northing values
easting <- as.numeric(coordinates[, 1])
northing <- as.numeric(coordinates[, 2])
# Convert UTM to lonlat using oce's utm2lonlat function
lonlat <- utm2lonlat(easting = easting, northing = northing, zone = 32, hemisphere = "N")
# Add the lonlat values as new columns
sequenced_individuals$longitude <- lonlat$longitude
sequenced_individuals$latitude <- lonlat$latitude
saveRDS(sequenced_individuals, file = "sequenced_individual_detailed_metadata.rds")



# Distance matrix for geographic locations ------
num_cores <- 8

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Create a function to calculate Haversine distance for a pair of coordinates
calculate_haversine_distance <- function(coord1, coord2) {
  distHaversine(coord1, coord2, r = 6378137)
}

# Extract latitude and longitude once
coordinates <- sampled_data[, c("longitude", "latitude")]

# Use matrix operations to parallelize the computation
geographic_distance_matrix_parallel <- foreach(i = 1:nrow(sampled_data), .combine = rbind) %dopar% {
  # Load the geosphere library inside the parallel block (only load once)
  library(geosphere)

  sapply(1:nrow(sampled_data), function(j) {
    calculate_haversine_distance(coordinates[i, ], coordinates[j, ])
  })
}

# Stop the parallel backend
stopCluster(cl)


# Now, geographic_distance_matrix_parallel contains the distance between individuals based on lat/long coordETRS89 coordinates using parallel processing
saveRDS(geographic_distance_matrix_parallel, file = "geographic_distance_matrix_subset.rds")
write.csv(geographic_distance_matrix_parallel, file = "geographic_distance_matrix_subset.csv")

# Matrix of differences in time between people -----
calculate_time_difference <- function(date1, date2) {
  as.numeric(difftime(date1, date2, units = "days"))
}
cl <- makeCluster(num_cores)
registerDoParallel(cl)
time_distance_matrix_parallel <- foreach(i = 1:nrow(sampled_data), .combine = rbind) %dopar% {
  sapply(1:nrow(sampled_data), function(j) {
    calculate_time_difference(sampled_data$date[i], sampled_data$date[j])
  })
}
stopCluster(cl)
# time_distance_matrix_df <- as.data.frame(time_distance_matrix_parallel)
saveRDS(time_distance_matrix_parallel, file = "time_distance_matrix.rds")
# write.csv(time_distance_matrix_parallel, file = "time_distance_matrix.csv")
