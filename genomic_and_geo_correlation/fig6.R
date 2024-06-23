# Geographic and Genomic Distance Correlation Plots

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr,
  lubridate, rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica,
  DSTdb, DBI, parallel, ggsignif, Rcpp, geosphere, biglm, graphics,
  pheatmap, viridis, patchwork, coefplot, ggpubr, gridExtra, hdf5r,
  JuliaCall
)


# Functions ------

# Function to extract upper triangular part of a matrix
lower_triangular <- function(mat) {
  mat[lower.tri(mat) | diag(rep(TRUE, nrow(mat)))] <- NA
  return(mat)
}

upper_triangular <- function(mat) {
  mat[upper.tri(mat) | diag(rep(TRUE, nrow(mat)))] <- NA
  return(mat)
}

read.jld2 <- function(fname, key) {
  file.h5 <- hdf5r::H5File$new(fname, mode = "r+")
  ds <- file.h5[[key]]
  dist <- array(ds[, ], dim = ds$dims)

  file.h5$close_all()

  return(dist)
}

# Load data --------

# Loading distance and time matrices from Open Street Maps for all datasets ----------------------
OSM_distances_full <- read.jld2()
OSM_distances_byzone <- read.jld2()
OSM_distances_landzone <- read.jld2()
OSM_distances_hovedstaden <- read.jld2()
OSM_distances_midtjylland <- read.jld2()
OSM_distances_nordjylland <- read.jld2()
OSM_distances_sjælland <- read.jld2()
OSM_distances_syddanmark <- read.jld2()
OSM_distances_copenhagen <- read.jld2()
OSM_distances_odense <- read.jld2()
OSM_distances_aarhus <- read.jld2()

OSM_distances_full_vector <- as.vector(lower_triangular(OSM_distances_full)) / 10000
OSM_distances_byzone_vector <- as.vector(lower_triangular(OSM_distances_byzone)) / 10000
OSM_distances_landzone_vector <- as.vector(lower_triangular(OSM_distances_landzone)) / 10000
OSM_distances_hovedstaden_vector <- as.vector(lower_triangular(OSM_distances_hovedstaden)) / 10000
OSM_distances_midtjylland_vector <- as.vector(lower_triangular(OSM_distances_midtjylland)) / 10000
OSM_distances_nordjylland_vector <- as.vector(lower_triangular(OSM_distances_nordjylland)) / 10000
OSM_distances_sjælland_vector <- as.vector(lower_triangular(OSM_distances_sjælland)) / 10000
OSM_distances_syddanmark_vector <- as.vector(lower_triangular(OSM_distances_syddanmark)) / 10000
OSM_distances_copenhagen_vector <- as.vector(lower_triangular(OSM_distances_copenhagen)) / 10000
OSM_distances_odense_vector <- as.vector(lower_triangular(OSM_distances_odense)) / 10000
OSM_distances_aarhus_vector <- as.vector(lower_triangular(OSM_distances_aarhus)) / 10000

# All of Denmark, 20,000 samples ---------
# Loading basic data
time_distance_matrix_parallel <- readRDS(file = "")
geographic_distance_matrix_parallel <- readRDS(file = "")
sequenced_individuals <- readRDS(file = "")
cophenetic_distances <- readRDS(file = "")
final_tree <- read.tree("")
final_distance_tree <- read.tree("")

cophenetic_distances_vector <- as.vector(lower_triangular(cophenetic_distances)) * 29891
time_distance_matrix_vector <- abs(as.vector(lower_triangular(time_distance_matrix_parallel))) / 7
geographic_distance_matrix_parallel_vector <- as.vector(lower_triangular(geographic_distance_matrix_parallel)) / 10000

combined_data <- data.frame(
  cophenetic_distances_vector = cophenetic_distances_vector,
  geographic_distance_matrix_parallel_vector = geographic_distance_matrix_parallel_vector,
  time_distance_matrix_vector = time_distance_matrix_vector,
  OSM_distances_full_vector = OSM_distances_full_vector
)

# Experiment ----

# Loading hamming distance matrix for national subsample
hamming_distances_national <- read.csv("", header = FALSE)
hamming_distances_national <- as.matrix(hamming_distances_national)
hamming_distances_national_vector <- as.vector(hamming_distances_national)

hamming_cophenetic_distances_vector <- as.vector(upper_triangular(cophenetic_distances)) * 29891
hamming_time_distance_matrix_vector <- abs(as.vector(upper_triangular(time_distance_matrix_parallel))) / 7
hamming_geographic_distance_matrix_parallel_vector <- as.vector(upper_triangular(geographic_distance_matrix_parallel)) / 10000
hamming_OSM_distances_full_vector <- as.vector(upper_triangular(OSM_distances_full)) / 10000

small_dataframe <- data.frame(
  hamming_cophenetic_distances_vector = hamming_cophenetic_distances_vector,
  hamming_geographic_distance_matrix_parallel_vector = hamming_geographic_distance_matrix_parallel_vector,
  hamming_distances_national_vector = hamming_distances_national_vector
)

hamming_geography_genetic_time_model <- biglm(hamming_distances_national_vector ~ hamming_geographic_distance_matrix_parallel_vector, data = small_dataframe)
summary(hamming_geography_genetic_time_model)

hamming_cophenetic_model <- biglm(hamming_distances_national_vector ~ hamming_cophenetic_distances_vector, data = small_dataframe)
summary(hamming_cophenetic_model)

# Regional data --------
# Regression by region, 10,000 samples for each
# Nordjylland
time_distance_matrix_Nordjylland <- readRDS(file = "")
geographic_distance_matrix_Nordjylland <- readRDS(file = "")
cophenetic_distances_Nordjylland <- readRDS(file = "")

geographic_distance_matrix_parallel_vector_Nordjylland <- as.vector(lower_triangular(geographic_distance_matrix_Nordjylland) / 10000)
time_distance_matrix_Nordjylland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Nordjylland)) / 7)
cophenetic_distances_vector_Nordjylland <- as.vector(lower_triangular(cophenetic_distances_Nordjylland) * 29891)

combined_data_Nordjylland <- data.frame(
  time_distance_matrix_Nordjylland_vector = time_distance_matrix_Nordjylland_vector,
  cophenetic_distances_vector_Nordjylland = cophenetic_distances_vector_Nordjylland,
  geographic_distance_matrix_parallel_vector_Nordjylland = geographic_distance_matrix_parallel_vector_Nordjylland,
  OSM_distances_nordjylland_vector = OSM_distances_nordjylland_vector
)

# Midtjylland
time_distance_matrix_Midtjylland <- readRDS(file = "")
geographic_distance_matrix_Midtjylland <- readRDS(file = "")
cophenetic_distances_Midtjylland <- readRDS(file = "")

geographic_distance_matrix_parallel_vector_Midtjylland <- as.vector(lower_triangular(geographic_distance_matrix_Midtjylland) / 10000)
time_distance_matrix_Midtjylland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Midtjylland)) / 7)
cophenetic_distances_vector_Midtjylland <- as.vector(lower_triangular(cophenetic_distances_Midtjylland) * 29891)

combined_data_Midtjylland <- data.frame(
  time_distance_matrix_Midtjylland_vector = time_distance_matrix_Midtjylland_vector,
  cophenetic_distances_vector_Midtjylland = cophenetic_distances_vector_Midtjylland,
  geographic_distance_matrix_parallel_vector_Midtjylland = geographic_distance_matrix_parallel_vector_Midtjylland,
  OSM_distances_midtjylland_vector = OSM_distances_midtjylland_vector
)

# Syddanmark
time_distance_matrix_Syddanmark <- readRDS(file = "")
geographic_distance_matrix_Syddanmark <- readRDS(file = "")
cophenetic_distances_Syddanmark <- readRDS(file = "")

geographic_distance_matrix_parallel_vector_Syddanmark <- as.vector(lower_triangular(geographic_distance_matrix_Syddanmark) / 10000)
time_distance_matrix_Syddanmark_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Syddanmark)) / 7)
cophenetic_distances_vector_Syddanmark <- as.vector(lower_triangular(cophenetic_distances_Syddanmark) * 29891)

combined_data_Syddanmark <- data.frame(
  time_distance_matrix_Syddanmark_vector = time_distance_matrix_Syddanmark_vector,
  cophenetic_distances_vector_Syddanmark = cophenetic_distances_vector_Syddanmark,
  geographic_distance_matrix_parallel_vector_Syddanmark = geographic_distance_matrix_parallel_vector_Syddanmark,
  OSM_distances_syddanmark_vector = OSM_distances_syddanmark_vector
)

# Hovedstaden
time_distance_matrix_Hovedstaden <- readRDS(file = "")
geographic_distance_matrix_Hovedstaden <- readRDS(file = "")
cophenetic_distances_Hovedstaden <- readRDS(file = "")

geographic_distance_matrix_parallel_vector_Hovedstaden <- as.vector(lower_triangular(geographic_distance_matrix_Hovedstaden) / 10000)
time_distance_matrix_Hovedstaden_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Hovedstaden)) / 7)
cophenetic_distances_vector_Hovedstaden <- as.vector(lower_triangular(cophenetic_distances_Hovedstaden) * 29891)

combined_data_Hovedstaden <- data.frame(
  time_distance_matrix_Hovedstaden_vector = time_distance_matrix_Hovedstaden_vector,
  cophenetic_distances_vector_Hovedstaden = cophenetic_distances_vector_Hovedstaden,
  geographic_distance_matrix_parallel_vector_Hovedstaden = geographic_distance_matrix_parallel_vector_Hovedstaden,
  OSM_distances_hovedstaden_vector = OSM_distances_hovedstaden_vector
)

# Sjælland
time_distance_matrix_Sjælland <- readRDS(file = "")
geographic_distance_matrix_Sjælland <- readRDS(file = "")
cophenetic_distances_Sjælland <- readRDS(file = "")

geographic_distance_matrix_parallel_vector_Sjælland <- as.vector(lower_triangular(geographic_distance_matrix_Sjælland) / 10000)
time_distance_matrix_Sjælland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Sjælland)) / 7)
cophenetic_distances_vector_Sjælland <- as.vector(lower_triangular(cophenetic_distances_Sjælland) * 29891)

combined_data_Sjælland <- data.frame(
  time_distance_matrix_Sjælland_vector = time_distance_matrix_Sjælland_vector,
  cophenetic_distances_vector_Sjælland = cophenetic_distances_vector_Sjælland,
  geographic_distance_matrix_parallel_vector_Sjælland = geographic_distance_matrix_parallel_vector_Sjælland,
  OSM_distances_sjælland_vector = OSM_distances_sjælland_vector
)



# Zone data-----------
# Byzone / Urban
sampled_data_reordered <- readRDS(file = "")
cophenetic_distances_byzone <- readRDS(file = "")
geographic_distance_matrix_byzone <- readRDS(file = "")
time_distance_matrix_byzone <- readRDS(file = "")

cophenetic_distances_byzone_vector <- as.vector(lower_triangular(cophenetic_distances_byzone)) * 29891
time_distance_matrix_byzone_vector <- abs(as.vector(lower_triangular(time_distance_matrix_byzone))) / 7
geographic_distance_matrix_byzone_vector <- as.vector(lower_triangular(geographic_distance_matrix_byzone)) / 10000

combined_data_byzone <- data.frame(
  cophenetic_distances_byzone_vector = cophenetic_distances_byzone_vector,
  geographic_distance_matrix_byzone_vector = geographic_distance_matrix_byzone_vector,
  time_distance_matrix_byzone_vector = time_distance_matrix_byzone_vector,
  OSM_distances_byzone_vector = OSM_distances_byzone_vector
)


# Landzone / Countryside
sampled_data_reordered_landzone <- readRDS(file = "")
cophenetic_distances_landzone <- readRDS(file = "")
geographic_distance_matrix_landzone <- readRDS(file = "")
time_distance_matrix_landzone <- readRDS(file = "")

cophenetic_distances_landzone_vector <- as.vector(lower_triangular(cophenetic_distances_landzone)) * 29891
time_distance_matrix_landzone_vector <- abs(as.vector(lower_triangular(time_distance_matrix_landzone))) / 7
geographic_distance_matrix_landzone_vector <- as.vector(lower_triangular(geographic_distance_matrix_landzone)) / 10000

combined_data_landzone <- data.frame(
  cophenetic_distances_landzone_vector = cophenetic_distances_landzone_vector,
  geographic_distance_matrix_landzone_vector = geographic_distance_matrix_landzone_vector,
  time_distance_matrix_landzone_vector = time_distance_matrix_landzone_vector,
  OSM_distances_landzone_vector = OSM_distances_landzone_vector
)


# Copenhagen
sampled_data_reordered_copenhagen <- readRDS(file = "")
cophenetic_distances_copenhagen <- readRDS(file = "")
geographic_distance_matrix_copenhagen <- readRDS(file = "")
time_distance_matrix_copenhagen <- readRDS(file = "")

cophenetic_distances_copenhagen_vector <- as.vector(lower_triangular(cophenetic_distances_copenhagen)) * 29891
time_distance_matrix_copenhagen_vector <- abs(as.vector(lower_triangular(time_distance_matrix_copenhagen))) / 7
geographic_distance_matrix_copenhagen_vector <- as.vector(lower_triangular(geographic_distance_matrix_copenhagen)) / 10000

combined_data_copenhagen <- data.frame(
  cophenetic_distances_copenhagen_vector = cophenetic_distances_copenhagen_vector,
  geographic_distance_matrix_copenhagen_vector = geographic_distance_matrix_copenhagen_vector,
  time_distance_matrix_copenhagen_vector = time_distance_matrix_copenhagen_vector,
  OSM_distances_copenhagen_vector = OSM_distances_copenhagen_vector
)

# Aarhus
sampled_data_reordered_aarhus <- readRDS(file = "")
cophenetic_distances_aarhus <- readRDS(file = "")
geographic_distance_matrix_aarhus <- readRDS(file = "")
time_distance_matrix_aarhus <- readRDS(file = "")

cophenetic_distances_aarhus_vector <- as.vector(lower_triangular(cophenetic_distances_aarhus)) * 29891
time_distance_matrix_aarhus_vector <- abs(as.vector(lower_triangular(time_distance_matrix_aarhus))) / 7
geographic_distance_matrix_aarhus_vector <- as.vector(lower_triangular(geographic_distance_matrix_aarhus)) / 10000

combined_data_aarhus <- data.frame(
  cophenetic_distances_aarhus_vector = cophenetic_distances_aarhus_vector,
  geographic_distance_matrix_aarhus_vector = geographic_distance_matrix_aarhus_vector,
  time_distance_matrix_aarhus_vector = time_distance_matrix_aarhus_vector,
  OSM_distances_aarhus_vector = OSM_distances_aarhus_vector
)


# Odense
sampled_data_reordered_odense <- readRDS(file = "")
cophenetic_distances_odense <- readRDS(file = "")
geographic_distance_matrix_odense <- readRDS(file = "")
time_distance_matrix_odense <- readRDS(file = "")

cophenetic_distances_odense_vector <- as.vector(lower_triangular(cophenetic_distances_odense)) * 29891
time_distance_matrix_odense_vector <- abs(as.vector(lower_triangular(time_distance_matrix_odense))) / 7
geographic_distance_matrix_odense_vector <- as.vector(lower_triangular(geographic_distance_matrix_odense)) / 10000

combined_data_odense <- data.frame(
  cophenetic_distances_odense_vector = cophenetic_distances_odense_vector,
  geographic_distance_matrix_odense_vector = geographic_distance_matrix_odense_vector,
  time_distance_matrix_odense_vector = time_distance_matrix_odense_vector,
  OSM_distances_odense_vector = OSM_distances_odense_vector
)






# Linear regression models, with time -------------

# Full sample
geography_genetic_time_model <- biglm(cophenetic_distances_vector ~ geographic_distance_matrix_parallel_vector + time_distance_matrix_vector, data = combined_data)
summary(geography_genetic_time_model)

OSM_genetic_time_model <- biglm(cophenetic_distances_vector ~ OSM_distances_full_vector + time_distance_matrix_vector, data = combined_data)
summary(OSM_genetic_time_model)

# Hovedstaden
geography_genetic_time_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ geographic_distance_matrix_parallel_vector_Hovedstaden + time_distance_matrix_Hovedstaden_vector, data = combined_data_Hovedstaden)
summary(geography_genetic_time_model_Hovedstaden)

OSM_genetic_time_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ OSM_distances_hovedstaden_vector + time_distance_matrix_Hovedstaden_vector, data = combined_data_Hovedstaden)
summary(OSM_genetic_time_model_Hovedstaden)

# Midtjylland
geography_genetic_time_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ geographic_distance_matrix_parallel_vector_Midtjylland + time_distance_matrix_Midtjylland_vector, data = combined_data_Midtjylland)
summary(geography_genetic_time_model_Midtjylland)

OSM_genetic_time_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ OSM_distances_midtjylland_vector + time_distance_matrix_Midtjylland_vector, data = combined_data_Midtjylland)
summary(OSM_genetic_time_model_Midtjylland)

# Nordjylland
geography_genetic_time_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ geographic_distance_matrix_parallel_vector_Nordjylland + time_distance_matrix_Nordjylland_vector, data = combined_data_Nordjylland)
summary(geography_genetic_time_model_Nordjylland)

OSM_genetic_time_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ OSM_distances_nordjylland_vector + time_distance_matrix_Nordjylland_vector, data = combined_data_Nordjylland)
summary(OSM_genetic_time_model_Nordjylland)

# Sjælland
geography_genetic_time_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ geographic_distance_matrix_parallel_vector_Sjælland + time_distance_matrix_Sjælland_vector, data = combined_data_Sjælland)
summary(geography_genetic_time_model_Sjælland)

OSM_genetic_time_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ OSM_distances_sjælland_vector + time_distance_matrix_Sjælland_vector, data = combined_data_Sjælland)
summary(OSM_genetic_time_model_Sjælland)

# Syddanmark
geography_genetic_time_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ geographic_distance_matrix_parallel_vector_Syddanmark + time_distance_matrix_Syddanmark_vector, data = combined_data_Syddanmark)
summary(geography_genetic_time_model_Syddanmark)

OSM_genetic_time_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ OSM_distances_syddanmark_vector + time_distance_matrix_Syddanmark_vector, data = combined_data_Syddanmark)
summary(OSM_genetic_time_model_Syddanmark)

# Byzone / Urban
geography_genetic_time_model_byzone <- biglm(cophenetic_distances_byzone_vector ~ geographic_distance_matrix_byzone_vector + time_distance_matrix_byzone_vector, data = combined_data_byzone)
summary(geography_genetic_time_model_byzone)

OSM_genetic_time_model_byzone <- biglm(cophenetic_distances_byzone_vector ~ OSM_distances_byzone_vector + time_distance_matrix_byzone_vector, data = combined_data_byzone)
summary(OSM_genetic_time_model_byzone)

# Landzone
geography_genetic_time_model_landzone <- biglm(cophenetic_distances_landzone_vector ~ geographic_distance_matrix_landzone_vector + time_distance_matrix_landzone_vector, data = combined_data_landzone)
summary(geography_genetic_time_model_landzone)

OSM_genetic_time_model_landzone <- biglm(cophenetic_distances_landzone_vector ~ OSM_distances_landzone_vector + time_distance_matrix_landzone_vector, data = combined_data_landzone)
summary(OSM_genetic_time_model_landzone)

# Copenhagen
geography_genetic_time_model_copenhagen <- biglm(cophenetic_distances_copenhagen_vector ~ geographic_distance_matrix_copenhagen_vector + time_distance_matrix_copenhagen_vector, data = combined_data_copenhagen)
summary(geography_genetic_time_model_copenhagen)

OSM_genetic_time_model_copenhagen <- biglm(cophenetic_distances_copenhagen_vector ~ OSM_distances_copenhagen_vector + time_distance_matrix_copenhagen_vector, data = combined_data_copenhagen)
summary(OSM_genetic_time_model_copenhagen)




# Linear regression models, without time -------------

# Full sample
geography_genetic_model <- biglm(cophenetic_distances_vector ~ geographic_distance_matrix_parallel_vector, data = combined_data)
summary(geography_genetic_model)

OSM_genetic_model <- biglm(cophenetic_distances_vector ~ OSM_distances_full_vector, data = combined_data)
summary(OSM_genetic_model)

# Hovedstaden
geography_genetic_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ geographic_distance_matrix_parallel_vector_Hovedstaden, data = combined_data_Hovedstaden)
summary(geography_genetic_model_Hovedstaden)

OSM_genetic_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ OSM_distances_hovedstaden_vector, data = combined_data_Hovedstaden)
summary(OSM_genetic_model_Hovedstaden)

# Midtjylland
geography_genetic_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ geographic_distance_matrix_parallel_vector_Midtjylland, data = combined_data_Midtjylland)
summary(geography_genetic_model_Midtjylland)

OSM_genetic_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ OSM_distances_midtjylland_vector, data = combined_data_Midtjylland)
summary(OSM_genetic_model_Midtjylland)

# Nordjylland
geography_genetic_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ geographic_distance_matrix_parallel_vector_Nordjylland, data = combined_data_Nordjylland)
summary(geography_genetic_model_Nordjylland)

OSM_genetic_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ OSM_distances_nordjylland_vector, data = combined_data_Nordjylland)
summary(OSM_genetic_model_Nordjylland)

# Sjælland
geography_genetic_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ geographic_distance_matrix_parallel_vector_Sjælland, data = combined_data_Sjælland)
summary(geography_genetic_model_Sjælland)

OSM_genetic_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ OSM_distances_sjælland_vector, data = combined_data_Sjælland)
summary(OSM_genetic_model_Sjælland)

# Syddanmark
geography_genetic_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ geographic_distance_matrix_parallel_vector_Syddanmark, data = combined_data_Syddanmark)
summary(geography_genetic_model_Syddanmark)

OSM_genetic_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ OSM_distances_syddanmark_vector, data = combined_data_Syddanmark)
summary(OSM_genetic_model_Syddanmark)

# Byzone / Urban
geography_genetic_model_byzone <- biglm(cophenetic_distances_byzone_vector ~ geographic_distance_matrix_byzone_vector, data = combined_data_byzone)
summary(geography_genetic_model_byzone)

OSM_genetic_model_byzone <- biglm(cophenetic_distances_byzone_vector ~ OSM_distances_byzone_vector, data = combined_data_byzone)
summary(OSM_genetic_model_byzone)

# Landzone
geography_genetic_model_landzone <- biglm(cophenetic_distances_landzone_vector ~ geographic_distance_matrix_landzone_vector, data = combined_data_landzone)
summary(geography_genetic_model_landzone)

OSM_genetic_model_landzone <- biglm(cophenetic_distances_landzone_vector ~ OSM_distances_landzone_vector, data = combined_data_landzone)
summary(OSM_genetic_model_landzone)

# Copenhagen
geography_genetic_model_copenhagen <- biglm(cophenetic_distances_copenhagen_vector ~ geographic_distance_matrix_copenhagen_vector, data = combined_data_copenhagen)
summary(geography_genetic_model_copenhagen)

OSM_genetic_model_copenhagen <- biglm(cophenetic_distances_copenhagen_vector ~ OSM_distances_copenhagen_vector, data = combined_data_copenhagen)
summary(OSM_genetic_model_copenhagen)






# Extracting coefficient values, with time -----------

# Full sample
coefficients <- coef(geography_genetic_time_model)
standard_errors <- sqrt(diag(vcov(geography_genetic_time_model)))

OSM_full_coefficients <- coef(OSM_genetic_time_model)
OSM_full_standard_errors <- sqrt(diag(vcov(OSM_genetic_time_model)))

# Hovedstaden
coefficients_Hovedstaden <- coef(geography_genetic_time_model_Hovedstaden)
standard_errors_Hovedstaden <- sqrt(diag(vcov(geography_genetic_time_model_Hovedstaden)))

OSM_coefficients_Hovedstaden <- coef(OSM_genetic_time_model_Hovedstaden)
OSM_standard_errors_Hovedstaden <- sqrt(diag(vcov(OSM_genetic_time_model_Hovedstaden)))

# Midtjylland
coefficients_Midtjylland <- coef(geography_genetic_time_model_Midtjylland)
standard_errors_Midtjylland <- sqrt(diag(vcov(geography_genetic_time_model_Midtjylland)))

OSM_coefficients_Midtjylland <- coef(OSM_genetic_time_model_Midtjylland)
OSM_standard_errors_Midtjylland <- sqrt(diag(vcov(OSM_genetic_time_model_Midtjylland)))

# Nordjylland
coefficients_Nordjylland <- coef(geography_genetic_time_model_Nordjylland)
standard_errors_Nordjylland <- sqrt(diag(vcov(geography_genetic_time_model_Nordjylland)))

OSM_coefficients_Nordjylland <- coef(OSM_genetic_time_model_Nordjylland)
OSM_standard_errors_Nordjylland <- sqrt(diag(vcov(OSM_genetic_time_model_Nordjylland)))

# Sjælland
coefficients_Sjælland <- coef(geography_genetic_time_model_Sjælland)
standard_errors_Sjælland <- sqrt(diag(vcov(geography_genetic_time_model_Sjælland)))

OSM_coefficients_Sjælland <- coef(OSM_genetic_time_model_Sjælland)
OSM_standard_errors_Sjælland <- sqrt(diag(vcov(OSM_genetic_time_model_Sjælland)))

# Syddanmark
coefficients_Syddanmark <- coef(geography_genetic_time_model_Syddanmark)
standard_errors_Syddanmark <- sqrt(diag(vcov(geography_genetic_time_model_Syddanmark)))

OSM_coefficients_Syddanmark <- coef(OSM_genetic_time_model_Syddanmark)
OSM_standard_errors_Syddanmark <- sqrt(diag(vcov(OSM_genetic_time_model_Syddanmark)))

# Byzone / Urban
coefficients_byzone <- coef(geography_genetic_time_model_byzone)
standard_errors_byzone <- sqrt(diag(vcov(geography_genetic_time_model_byzone)))

OSM_coefficients_byzone <- coef(OSM_genetic_time_model_byzone)
OSM_standard_errors_byzone <- sqrt(diag(vcov(OSM_genetic_time_model_byzone)))

# Landzone
coefficients_landzone <- coef(geography_genetic_time_model_landzone)
standard_errors_landzone <- sqrt(diag(vcov(geography_genetic_time_model_landzone)))

OSM_coefficients_landzone <- coef(OSM_genetic_time_model_landzone)
OSM_standard_errors_landzone <- sqrt(diag(vcov(OSM_genetic_time_model_landzone)))

# Copenhagen
coefficients_copenhagen <- coef(geography_genetic_time_model_copenhagen)
standard_errors_copenhagen <- sqrt(diag(vcov(geography_genetic_time_model_copenhagen)))

OSM_coefficients_copenhagen <- coef(OSM_genetic_time_model_copenhagen)
OSM_standard_errors_copenhagen <- sqrt(diag(vcov(OSM_genetic_time_model_copenhagen)))



# Extracting coefficient values, without time -----------

# Full sample
coefficients_no_time <- coef(geography_genetic_model)
standard_errors_no_time <- sqrt(diag(vcov(geography_genetic_model)))

OSM_full_coefficients_no_time <- coef(OSM_genetic_model)
OSM_full_standard_errors_no_time <- sqrt(diag(vcov(OSM_genetic_model)))

# Hovedstaden
coefficients_Hovedstaden_no_time <- coef(geography_genetic_model_Hovedstaden)
standard_errors_Hovedstaden_no_time <- sqrt(diag(vcov(geography_genetic_model_Hovedstaden)))

OSM_coefficients_Hovedstaden_no_time <- coef(OSM_genetic_model_Hovedstaden)
OSM_standard_errors_Hovedstaden_no_time <- sqrt(diag(vcov(OSM_genetic_model_Hovedstaden)))

# Midtjylland
coefficients_Midtjylland_no_time <- coef(geography_genetic_model_Midtjylland)
standard_errors_Midtjylland_no_time <- sqrt(diag(vcov(geography_genetic_model_Midtjylland)))

OSM_coefficients_Midtjylland_no_time <- coef(OSM_genetic_model_Midtjylland)
OSM_standard_errors_Midtjylland_no_time <- sqrt(diag(vcov(OSM_genetic_model_Midtjylland)))

# Nordjylland
coefficients_Nordjylland_no_time <- coef(geography_genetic_model_Nordjylland)
standard_errors_Nordjylland_no_time <- sqrt(diag(vcov(geography_genetic_model_Nordjylland)))

OSM_coefficients_Nordjylland_no_time <- coef(OSM_genetic_model_Nordjylland)
OSM_standard_errors_Nordjylland_no_time <- sqrt(diag(vcov(OSM_genetic_model_Nordjylland)))

# Sjælland
coefficients_Sjælland_no_time <- coef(geography_genetic_model_Sjælland)
standard_errors_Sjælland_no_time <- sqrt(diag(vcov(geography_genetic_model_Sjælland)))

OSM_coefficients_Sjælland_no_time <- coef(OSM_genetic_model_Sjælland)
OSM_standard_errors_Sjælland_no_time <- sqrt(diag(vcov(OSM_genetic_model_Sjælland)))

# Syddanmark
coefficients_Syddanmark_no_time <- coef(geography_genetic_model_Syddanmark)
standard_errors_Syddanmark_no_time <- sqrt(diag(vcov(geography_genetic_model_Syddanmark)))

OSM_coefficients_Syddanmark_no_time <- coef(OSM_genetic_model_Syddanmark)
OSM_standard_errors_Syddanmark_no_time <- sqrt(diag(vcov(OSM_genetic_model_Syddanmark)))

# Byzone / Urban
coefficients_byzone_no_time <- coef(geography_genetic_model_byzone)
standard_errors_byzone_no_time <- sqrt(diag(vcov(geography_genetic_model_byzone)))

OSM_coefficients_byzone_no_time <- coef(OSM_genetic_model_byzone)
OSM_standard_errors_byzone_no_time <- sqrt(diag(vcov(OSM_genetic_model_byzone)))

# Landzone
coefficients_landzone_no_time <- coef(geography_genetic_model_landzone)
standard_errors_landzone_no_time <- sqrt(diag(vcov(geography_genetic_model_landzone)))

OSM_coefficients_landzone_no_time <- coef(OSM_genetic_model_landzone)
OSM_standard_errors_landzone_no_time <- sqrt(diag(vcov(OSM_genetic_model_landzone)))

# Copenhagen
coefficients_copenhagen_no_time <- coef(geography_genetic_model_copenhagen)
standard_errors_copenhagen_no_time <- sqrt(diag(vcov(geography_genetic_model_copenhagen)))

OSM_coefficients_copenhagen_no_time <- coef(OSM_genetic_model_copenhagen)
OSM_standard_errors_copenhagen_no_time <- sqrt(diag(vcov(OSM_genetic_model_copenhagen)))



# Combining coefficients, time -------

# Combine coefficients and standard errors into a data frame
geography_coefficients <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    coefficients[2], coefficients_Hovedstaden[2],
    coefficients_Midtjylland[2], coefficients_Nordjylland[2],
    coefficients_Sjælland[2], coefficients_Syddanmark[2],
    coefficients_byzone[2], coefficients_landzone[2], coefficients_copenhagen[2]
  ),
  std_error = c(
    standard_errors[2], standard_errors_Hovedstaden[2],
    standard_errors_Midtjylland[2], standard_errors_Nordjylland[2],
    standard_errors_Sjælland[2], standard_errors_Syddanmark[2],
    standard_errors_byzone[2], standard_errors_landzone[2], standard_errors_copenhagen[2]
  )
)
geography_coefficients$conf_low <- geography_coefficients$estimate - 1.96 * geography_coefficients$std_error
geography_coefficients$conf_high <- geography_coefficients$estimate + 1.96 * geography_coefficients$std_error



# OSM
OSM_coefficients <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    OSM_full_coefficients[2], OSM_coefficients_Hovedstaden[2],
    OSM_coefficients_Midtjylland[2], OSM_coefficients_Nordjylland[2],
    OSM_coefficients_Sjælland[2], OSM_coefficients_Syddanmark[2],
    OSM_coefficients_byzone[2], OSM_coefficients_landzone[2], OSM_coefficients_copenhagen[2]
  ),
  std_error = c(
    OSM_full_standard_errors[2], OSM_standard_errors_Hovedstaden[2],
    OSM_standard_errors_Midtjylland[2], OSM_standard_errors_Nordjylland[2],
    OSM_standard_errors_Sjælland[2], OSM_standard_errors_Syddanmark[2],
    OSM_standard_errors_byzone[2], OSM_standard_errors_landzone[2], OSM_standard_errors_copenhagen[2]
  )
)
OSM_coefficients$conf_low <- OSM_coefficients$estimate - 1.96 * OSM_coefficients$std_error
OSM_coefficients$conf_high <- OSM_coefficients$estimate + 1.96 * OSM_coefficients$std_error


# Time coefficients
geography_coefficients_time <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    coefficients[3], coefficients_Hovedstaden[3],
    coefficients_Midtjylland[3], coefficients_Nordjylland[3],
    coefficients_Sjælland[3], coefficients_Syddanmark[3],
    coefficients_byzone[3], coefficients_landzone[3], coefficients_copenhagen[3]
  ),
  std_error = c(
    standard_errors[3], standard_errors_Hovedstaden[3],
    standard_errors_Midtjylland[3], standard_errors_Nordjylland[3],
    standard_errors_Sjælland[3], standard_errors_Syddanmark[3],
    standard_errors_byzone[3], standard_errors_landzone[3], standard_errors_copenhagen[3]
  )
)
geography_coefficients_time$conf_low <- geography_coefficients_time$estimate - 1.96 * geography_coefficients_time$std_error
geography_coefficients_time$conf_high <- geography_coefficients_time$estimate + 1.96 * geography_coefficients_time$std_error


# For OSM Coefficients with only the third term
OSM_coefficients_time <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    OSM_full_coefficients[3], OSM_coefficients_Hovedstaden[3],
    OSM_coefficients_Midtjylland[3], OSM_coefficients_Nordjylland[3],
    OSM_coefficients_Sjælland[3], OSM_coefficients_Syddanmark[3],
    OSM_coefficients_byzone[3], OSM_coefficients_landzone[3], OSM_coefficients_copenhagen[3]
  ),
  std_error = c(
    OSM_full_standard_errors[3], OSM_standard_errors_Hovedstaden[3],
    OSM_standard_errors_Midtjylland[3], OSM_standard_errors_Nordjylland[3],
    OSM_standard_errors_Sjælland[3], OSM_standard_errors_Syddanmark[3],
    OSM_standard_errors_byzone[3], OSM_standard_errors_landzone[3], OSM_standard_errors_copenhagen[3]
  )
)
OSM_coefficients_time$conf_low <- OSM_coefficients_time$estimate - 1.96 * OSM_coefficients_time$std_error
OSM_coefficients_time$conf_high <- OSM_coefficients_time$estimate + 1.96 * OSM_coefficients_time$std_error







# Combining coefficients, without time -------

# Combine coefficients and standard errors into a data frame
geography_coefficients_no_time <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    coefficients_no_time[2], coefficients_Hovedstaden_no_time[2],
    coefficients_Midtjylland_no_time[2], coefficients_Nordjylland_no_time[2],
    coefficients_Sjælland_no_time[2], coefficients_Syddanmark_no_time[2],
    coefficients_byzone_no_time[2], coefficients_landzone_no_time[2], coefficients_copenhagen_no_time[2]
  ),
  std_error = c(
    standard_errors_no_time[2], standard_errors_Hovedstaden_no_time[2],
    standard_errors_Midtjylland_no_time[2], standard_errors_Nordjylland_no_time[2],
    standard_errors_Sjælland_no_time[2], standard_errors_Syddanmark_no_time[2],
    standard_errors_byzone_no_time[2], standard_errors_landzone_no_time[2], standard_errors_copenhagen_no_time[2]
  )
)
geography_coefficients_no_time$conf_low <- geography_coefficients_no_time$estimate - 1.96 * geography_coefficients_no_time$std_error
geography_coefficients_no_time$conf_high <- geography_coefficients_no_time$estimate + 1.96 * geography_coefficients_no_time$std_error



# OSM
OSM_coefficients_no_time <- data.frame(
  model = c(
    "National", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark",
    "Urban", "Countryside", "Copenhagen"
  ),
  estimate = c(
    OSM_full_coefficients_no_time[2], OSM_coefficients_Hovedstaden_no_time[2],
    OSM_coefficients_Midtjylland_no_time[2], OSM_coefficients_Nordjylland_no_time[2],
    OSM_coefficients_Sjælland_no_time[2], OSM_coefficients_Syddanmark_no_time[2],
    OSM_coefficients_byzone_no_time[2], OSM_coefficients_landzone_no_time[2], OSM_coefficients_copenhagen_no_time[2]
  ),
  std_error = c(
    OSM_full_standard_errors_no_time[2], OSM_standard_errors_Hovedstaden_no_time[2],
    OSM_standard_errors_Midtjylland_no_time[2], OSM_standard_errors_Nordjylland_no_time[2],
    OSM_standard_errors_Sjælland_no_time[2], OSM_standard_errors_Syddanmark_no_time[2],
    OSM_standard_errors_byzone_no_time[2], OSM_standard_errors_landzone_no_time[2], OSM_standard_errors_copenhagen_no_time[2]
  )
)
OSM_coefficients_no_time$conf_low <- OSM_coefficients_no_time$estimate - 1.96 * OSM_coefficients_no_time$std_error
OSM_coefficients_no_time$conf_high <- OSM_coefficients_no_time$estimate + 1.96 * OSM_coefficients_no_time$std_error















# Plotting coefficients, with time ---------
model_order <- c("National", "Urban", "Countryside", "Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark", "Copenhagen")
geography_coefficients$color <- ifelse(geography_coefficients$estimate > 0, "Above", "Below")

# For geographic plot
geography_coefficients$model <- factor(geography_coefficients$model, levels = model_order)
geographic_plot <- ggplot(geography_coefficients, aes(x = model, y = estimate, color = color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  scale_color_manual(values = c("Above" = "#440154FF", "Below" = "#21908CFF")) +
  labs(y = "Molecular Change", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  ggtitle("Per 10km, Euclidean") +
  guides(color = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

geographic_plot

# For geographic plot
OSM_coefficients$model <- factor(OSM_coefficients$model, levels = model_order)
OSM_coefficients$color <- ifelse(OSM_coefficients$estimate > 0, "Above", "Below")

OSM_plot <- ggplot(OSM_coefficients, aes(x = model, y = estimate, color = color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  scale_color_manual(values = c("Above" = "#440154FF", "Below" = "#21908CFF")) +
  labs(y = "Molecular Change", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  ggtitle("Per 10km, Travel") +
  guides(color = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plotting coefficients, without time ---------
geography_coefficients_no_time$color <- ifelse(geography_coefficients_no_time$estimate > 0, "Above", "Below")

# For geographic plot
geography_coefficients_no_time$model <- factor(geography_coefficients_no_time$model, levels = model_order)
geographic_plot_no_time <- ggplot(geography_coefficients_no_time, aes(x = model, y = estimate, color = color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  scale_color_manual(values = c("Above" = "#440154FF", "Below" = "#21908CFF")) +
  labs(y = "Molecular Change", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  ggtitle("Per 10km, Euclidean") +
  guides(color = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# For geographic plot
OSM_coefficients_no_time$model <- factor(OSM_coefficients_no_time$model, levels = model_order)
OSM_coefficients_no_time$color <- ifelse(OSM_coefficients_no_time$estimate > 0, "Above", "Below")

OSM_plot_no_time <- ggplot(OSM_coefficients_no_time, aes(x = model, y = estimate, color = color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  scale_color_manual(values = c("Above" = "#440154FF", "Below" = "#21908CFF")) +
  labs(y = "Molecular Change", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  ggtitle("Per 10km, Travel") +
  guides(color = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Cophenetic distance plot -------

results_by_region <- readRDS(file = "results_by_region.rds")
total_sample_data <- readRDS(file = "total_sample_data.rds")
total_sample_data$region_name <- "Full"
total_sample_data$Region <- 1
merged_data <- rbind(results_by_region, total_sample_data)

cophenetic_mean_distance_plot <- ggplot(merged_data, aes(x = date, y = mean_distance * 29891, color = region_name)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette
  labs(title = NULL, x = "Date", y = "Mean Cophenetic Distance", color = "Region") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Getting legend to include in final plot
legend_plot <- ggplot(merged_data, aes(x = date, y = mean_distance * 29891, color = region_name)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette
  labs(title = "Mean Cophenetic Distance over Time", x = "Date", y = "Mean Cophenetic Distance", color = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
legend_plot <- get_legend(legend_plot)







# Household analysis -----------
distance_matrix_average_three_people <- readRDS(file = "")
distance_matrix_time_three_people <- readRDS(file = "")
geographic_distance_matrix_households <- readRDS(file = "")
samples_households <- readRDS(file = "")

# Histograms
extract_diagonal_and_non_diagonal <- function(matrix) {
  diagonal <- diag(matrix)
  non_diagonal <- matrix[lower.tri(matrix, diag = FALSE)]
  non_diagonal <- non_diagonal[!is.na(non_diagonal)]
  return(list(diagonal = diagonal, non_diagonal = non_diagonal))
}
# Extract diagonal and non-diagonal elements for both matrices
distance_matrix_average_three_people_nucleotides <- distance_matrix_average_three_people * 29891
elements_average <- extract_diagonal_and_non_diagonal(distance_matrix_average_three_people_nucleotides)

histogram_diagonal_unadjusted <- ggplot(data = data.frame(x = elements_average$diagonal), aes(x = x)) +
  geom_histogram(fill = "#440154", color = "white", bins = 30) +
  labs(title = "Within Households", x = "Mean Pairwise Cophenetic Distance", y = "Frequency") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white"))

# Histogram for the non-diagonal elements in average matrix
histogram_non_diagonal_unadjusted <- ggplot(data = data.frame(x = elements_average$non_diagonal), aes(x = x)) +
  geom_histogram(fill = "#440154", color = "white", bins = 30) +
  labs(title = "Between Households", x = "Mean Pairwise Cophenetic Distance", y = "Frequency") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(0, 90000))

# Normalizing for time
time_normalized_cophenetic_distance <- (distance_matrix_average_three_people / distance_matrix_time_three_people) * 29891
time_normalized_average_distance <- extract_diagonal_and_non_diagonal(time_normalized_cophenetic_distance)

# Histogram for the diagonal elements in average matrix with time
histogram_diagonal_adjusted <- ggplot(data = data.frame(x = time_normalized_average_distance$diagonal), aes(x = x)) +
  geom_histogram(fill = "#440154", color = "white", bins = 30) +
  labs(title = "Within Households, Time-Normalized", x = "Mean Pairwise Cophenetic Distance", y = "Frequency") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white")) + # Change background to plain
  scale_x_continuous(limits = c(-0.1, 3), breaks = seq(0, 3, by = 0.5)) # Set x limits and tick marks

# Histogram for the non-diagonal elements in average matrix with time
histogram_non_diagonal_adjusted <- ggplot(data = data.frame(x = time_normalized_average_distance$non_diagonal), aes(x = x)) +
  geom_histogram(fill = "#440154", color = "white", bins = 30) +
  labs(title = "Between Households, Time-Normalized", x = "Mean Pairwise Cophenetic Distance", y = "Frequency") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white")) + # Change background to plain
  scale_x_continuous(limits = c(-0.1, 3), breaks = seq(0, 3, by = 0.5)) +
  scale_y_continuous(limits = c(0, 90000))


# Linear regression model between households
household_cophenetic_distance_matrix_vector <- as.vector(upper_triangular(distance_matrix_average_three_people)) * 29891
household_time_distance_matrix_vector <- as.vector(upper_triangular(distance_matrix_time_three_people)) / 7
household_geographic_distance_matrix_households_vector <- as.vector(upper_triangular(geographic_distance_matrix_households)) / 10000

household_combined_data <- data.frame(
  cophenetic_distance_matrix_vector = household_cophenetic_distance_matrix_vector,
  geographic_distance_matrix_households_vector = household_geographic_distance_matrix_households_vector,
  time_distance_matrix_vector = household_time_distance_matrix_vector
)

household_geography_genetic_time_model <- biglm(cophenetic_distance_matrix_vector ~ geographic_distance_matrix_households_vector + time_distance_matrix_vector, data = household_combined_data)
summary(household_geography_genetic_time_model)

coefficients_household <- coef(household_geography_genetic_time_model)
standard_errors_household <- sqrt(diag(vcov(household_geography_genetic_time_model)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_household <- data.frame(
  term = c("Per 10km, Euclidean", "Per 1 Week"),
  estimate = c(coefficients_household[2], coefficients_household[3]),
  std_error = c(standard_errors_household[2], standard_errors_household[3])
)
coefs_ci_manual_household$conf_low <- coefs_ci_manual_household$estimate - 1.96 * coefs_ci_manual_household$std_error
coefs_ci_manual_household$conf_high <- coefs_ci_manual_household$estimate + 1.96 * coefs_ci_manual_household$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_household <- ggplot(coefs_ci_manual_household, aes(x = term, y = estimate, color = term)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = "Molecular Change", x = NULL) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  guides(color = FALSE, shape = FALSE) +
  scale_color_manual(values = c("#440154FF", "#21908CFF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

coefficients_plot_household
ggsave("", coefficients_plot_household, width = 9, height = 6)




# Overall, combined plot --------
combined_geo_plot <- grid.arrange(histogram_diagonal_unadjusted,
  histogram_non_diagonal_unadjusted, geographic_plot, OSM_plot,
  histogram_diagonal_adjusted,
  histogram_non_diagonal_adjusted, cophenetic_mean_distance_plot,
  legend_plot,
  ncol = 4
)

ggsave("", combined_geo_plot, width = 15, height = 7)
