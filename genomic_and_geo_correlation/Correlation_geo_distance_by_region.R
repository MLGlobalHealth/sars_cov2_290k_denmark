# Data generation for correlation between geographic region and cophenetic distance by region

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Reading in packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, ggplotify,
  phangorn, heatmaply, RColorBrewer, graphics, purrr, future.apply, geosphere, patchwork, coefplot,
  adephylo, biglm, pheatmap
)

# Data Preparation (Tree and Metadata) ------------------------------
sequenced_individuals <- readRDS(file = "data/sequenced_individuals.RDS")
distance_tree <- read.tree("data/tree/distance_tree.tree")

# Change tip labels to PERSON_ID values
strain_person_id_mapping <- setNames(sequenced_individuals$PERSON_ID, sequenced_individuals$strain)
matched_person_ids <- strain_person_id_mapping[distance_tree$tip.label]
distance_tree$tip.label <- matched_person_ids

# Function to calculate Haversine distance for a pair of coordinates
haversine_distance <- function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2), r = 6378137)
}

# Function to calculate distance matrix based on longitude and latitude using parallel processing
calculate_distance_matrix_parallel <- function(data) {
  num_individuals <- nrow(data)
  library(geosphere)
  distances <- foreach(
    i = 1:num_individuals,
    .combine = rbind
  ) %dopar% {
    result <- numeric(num_individuals)
    for (j in 1:num_individuals) {
      result[j] <- haversine_distance(
        data[i, "longitude"],
        data[i, "latitude"],
        data[j, "longitude"],
        data[j, "latitude"]
      )
    }
    return(result)
  }

  return(distances)
}

# Define region names and corresponding REGIONSKODE values
region_names <- c("Region Nordjylland", "Region Midtjylland", "Region Syddanmark", "Region Hovedstaden", "Region Sjælland")
region_codes <- c(1081, 1082, 1083, 1084, 1085)

# Loop through each region
for (i in seq_along(region_names)) {
  region_name <- region_names[i]
  region_code <- region_codes[i]

  # Filter data for the current region
  region_data <- subset(sequenced_individuals, REGIONSKODE == region_code)

  # Cophenetic distance matrix for 10,000 random samples
  set.seed(123)
  sampled_data <- data.frame()
  sampled_person_ids <- character(0) # Empty character vector to store sampled PERSON_ID values
  # Loop until you have sampled 10,000 unique rows
  while (nrow(sampled_data) < 10000) {
    # Sample a new row
    new_row <- region_data[sample(nrow(region_data), 1), ]

    if (!(new_row$PERSON_ID %in% sampled_person_ids)) {
      sampled_data <- rbind(sampled_data, new_row)
      sampled_person_ids <- c(sampled_person_ids, new_row$PERSON_ID)
    }
  }
  matched_strains <- match(distance_tree$tip.label, sampled_person_ids)
  # Identify tips to keep based on matched_strains
  tips_to_keep <- distance_tree$tip.label[!is.na(matched_strains)]
  tips_to_keep <- as.character(tips_to_keep)
  pruned_tree <- keep.tip(distance_tree, tips_to_keep)
  cophenetic_distances <- cophenetic(pruned_tree)
  # Save cophenetic distances for the region
  saveRDS(cophenetic_distances, file = paste0("cophenetic_distances_", gsub(" ", "_", region_name), ".rds"))

  order_of_names <- rownames(cophenetic_distances)
  sampled_data <- sampled_data[match(order_of_names, sampled_data$PERSON_ID), ]

  library(geosphere)
  # Convert UTM coordinates to longitude and latitude
  sampled_data$longitude_and_latitude <- paste(sampled_data$ETRS89_EAST.x, sampled_data$ETRS89_NORTH.x, sep = ", ")
  coordinates <- strsplit(as.character(sampled_data$longitude_and_latitude), ", ")
  coordinates <- matrix(unlist(coordinates), ncol = 2, byrow = TRUE)
  easting <- as.numeric(coordinates[, 1])
  northing <- as.numeric(coordinates[, 2])
  lonlat <- utm2lonlat(easting = easting, northing = northing, zone = 32, hemisphere = "N")
  sampled_data$longitude <- lonlat$longitude
  sampled_data$latitude <- lonlat$latitude
  saveRDS(sampled_data, file = paste0("sampled_data_", gsub(" ", "_", region_name), ".rds"))

  # Distance matrix for geographic locations using parallel processing
  num_cores <- 8
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  library(geosphere)
  # Export the haversine_distance function to the parallel workers
  clusterExport(cl, "haversine_distance")
  clusterExport(cl, "distHaversine")


  geographic_distance_matrix_parallel <- calculate_distance_matrix_parallel(sampled_data)

  saveRDS(geographic_distance_matrix_parallel, file = paste0("geographic_distance_matrix_", gsub(" ", "_", region_name), ".rds"))

  stopCluster(cl)
}


# Getting time matrices ----------
# Define the calculate_time_difference function
calculate_time_difference <- function(date1, date2) {
  as.numeric(difftime(date1, date2, units = "days"))
}

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Loop through each region
for (i in seq_along(region_names)) {
  region_name <- region_names[i]

  # Load sampled_data for the current region
  sampled_data <- readRDS(paste0("sampled_data_", gsub(" ", "_", region_name), ".rds"))
  sampled_data$date <- as.Date(sampled_data$date) # Ensure 'date' is in Date format

  time_distance_matrix_parallel <- foreach(i = seq_len(nrow(sampled_data)), .combine = rbind) %dopar% {
    sapply(seq_len(nrow(sampled_data)), function(j) {
      calculate_time_difference(sampled_data$date[i], sampled_data$date[j])
    })
  }

  saveRDS(time_distance_matrix_parallel, file = paste0("time_distance_matrix_", gsub(" ", "_", region_name), ".rds"))
}
stopCluster(cl)
