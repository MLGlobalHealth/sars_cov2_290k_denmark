# Castor analysis

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages ----------
library(
  ape, phytools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel,
  parallel, ggsignif, viridis, ggtree, gridExtra, cowplot, ggplotify,
  RColorBrewer, graphics, patchwork, coefplot, castor
)

# Loading basic data -----
sequenced_individuals <- readRDS(file = "")
final_tree <- read.tree("")
final_distance_tree <- read.tree("")
geo_distance_matrix_parallel <- readRDS(file="")
time_distance_matrix_parallel <- readRDS(file="")
cophenetic_distances <- readRDS(file="")

rownames(geo_distance_matrix_parallel) <- colnames(geo_distance_matrix_parallel) <- sequenced_individuals$PERSON_ID

rownames(time_distance_matrix_parallel) <- colnames(time_distance_matrix_parallel) <- sequenced_individuals$PERSON_ID

rownames(cophenetic_distances) <- colnames(cophenetic_distances) <- sequenced_individuals$PERSON_ID

metadata_sequenced_individuals <- readRDS(file = "")


# Correlation geo distance and phylogenetic distance using castor instead -----
# Getting latitudes and longitudes to input
tip_labels_order <- final_distance_tree$tip.label

filtered_data <- sequenced_individuals %>%
  filter(!is.na(latitude) & !is.na(longitude) & strain %in% tip_labels_order)

pruned_tree <- drop.tip(final_distance_tree, setdiff(tip_labels_order, filtered_data$strain))

new_tip_order <- pruned_tree$tip.label

tip_latitudes <- filtered_data$latitude[match(new_tip_order, filtered_data$strain)]

tip_longitudes <- filtered_data$longitude[match(new_tip_order, filtered_data$strain)]

phylo_geodistance_correlation <- correlate_phylo_geodistances(pruned_tree,
  tip_latitudes,
  tip_longitudes,
  correlation_method = "pearson",
  max_phylodistance = Inf,
  Npermutations = 1000,
  alternative = "two-sided",
  radius = 6371
)


phylo_geodistance_correlation <- readRDS(file = "")
# Results $correlation, [1] -0.00808918, $Npairs [1] 400000000, $Pvalue [1] 0, $mean_random_correlation [1] 0.0003436232



# Diffusion analysis ---------
tip_labels_order <- final_tree$tip.label
filtered_data <- sequenced_individuals %>%
  filter(!is.na(latitude) & !is.na(longitude) & strain %in% tip_labels_order)
pruned_time_tree <- drop.tip(final_tree, setdiff(tip_labels_order, filtered_data$strain))
new_time_tree_tip_order <- pruned_time_tree$tip.label
tip_time_tree_latitudes <- filtered_data$latitude[match(new_time_tree_tip_order, filtered_data$strain)]
tip_time_tree_longitudes <- filtered_data$longitude[match(new_time_tree_tip_order, filtered_data$strain)]

# Sampling fraction of 0.4
sampling_fraction <- 0.4
diffusion_results_sf <- castor::fit_sbm_geobiased_const(pruned_time_tree,
  tip_time_tree_latitudes,
  tip_time_tree_longitudes,
  radius = 6371,
  Nsims = 1000,
  Nbootstraps = 1000,
  rarefaction = sampling_fraction
)
saveRDS(
  diffusion_results_sf,
  file = paste0("diffusion_results_sf_", sampling_fraction, ".RDS", sep = ""),
)
