#Correlation plots
# Read relevant packages ----------
library(ape)
library(phytools)
library(TreeTools)
library(dplyr)
library(tidyverse)
library(data.table)
library(dbplyr)
library(lubridate)
library(rlang)
library(foreach)
library(doParallel)
library(DSTora)
library(ROracle)
library(DSTcolectica)
library(DSTdb)
library(DBI)
library(parallel)
library(doParallel)
library(foreach)
library(ggsignif)
library(Rcpp)
library(geosphere)
library(biglm)
library(graphics)
library(pheatmap)
library(viridis)
library(geosphere)
library(patchwork)
library(coefplot)
library(ggpubr)


# Functions ------
# Function to extract lower triangular part of a matrix
lower_triangular <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}


# All of Denmark, 20,000 samples  -------------

#Loading basic data -----
time_distance_matrix_parallel <- readRDS(file = "time_distance_matrix.rds")
geographic_distance_matrix_parallel <- readRDS(file = "geographic_distance_matrix_subset.rds")
sequenced_individuals <- readRDS(file="sampled_subset_data.rds")
cophenetic_distances <- readRDS(file="cophenetic_distances_subsample.rds")
final_tree <- read.tree("final_pruned_no_matOptimize_chronumental_consensus_2021_tree.tree")
final_distance_tree <- read.tree("pruned_newids_consensus_2021_tree.tree")
rownames(geographic_distance_matrix_parallel) <- colnames(geographic_distance_matrix_parallel) <- sequenced_individuals$PERSON_ID
rownames(time_distance_matrix_parallel) <- colnames(time_distance_matrix_parallel) <- sequenced_individuals$PERSON_ID
rownames(cophenetic_distances) <- colnames(cophenetic_distances) <- sequenced_individuals$PERSON_ID

# Apply the function to matrices
lower_tri_geographic <- lower_triangular(geographic_distance_matrix_parallel) / 10000
lower_tri_time <- abs(lower_triangular(time_distance_matrix_parallel)) / 7
lower_tri_cophenetic <- lower_triangular(cophenetic_distances) * 29891

# Flatten the lower triangular matrices to vectors
geographic_distance_matrix_parallel_vector <- as.vector(lower_tri_geographic)
time_distance_matrix_vector <- as.vector(lower_tri_time)
cophenetic_distances_vector <- as.vector(lower_tri_cophenetic)

combined_data <- data.frame(
  cophenetic_distances_vector = cophenetic_distances_vector,
  geographic_distance_matrix_parallel_vector = geographic_distance_matrix_parallel_vector,
  time_distance_matrix_vector = time_distance_matrix_vector)


# Regress genetic distance with geographic distance and time distance -----
geography_genetic_time_model <- biglm(cophenetic_distances_vector ~ geographic_distance_matrix_parallel_vector + time_distance_matrix_vector, data = combined_data)
summary(geography_genetic_time_model)

#With interaction term:
# Regress genetic distance with geographic distance, time distance, and their interaction term
geography_genetic_time_interaction_model <- biglm(cophenetic_distances_vector ~ geographic_distance_matrix_parallel_vector * time_distance_matrix_vector, data = combined_data)
summary(geography_genetic_time_interaction_model)

#Linear regression
# Get coefficients and their standard errors for the first model
coefficients_model1 <- coef(geography_genetic_time_model)
standard_errors_model1 <- sqrt(diag(vcov(geography_genetic_time_model)))

# Get coefficients and their standard errors for the second model
coefficients_model2 <- coef(geography_genetic_time_interaction_model)
standard_errors_model2 <- sqrt(diag(vcov(geography_genetic_time_interaction_model)))

# Combine coefficients and standard errors into a data frame
coefs_ci_manual <- data.frame(
  model = rep(c("Model 1", "Model 2"), each = 2),
  term = c("Geographic Distance", "Time Distance"),
  estimate = c(coefficients_model1[2], coefficients_model1[3],
               coefficients_model2[2], coefficients_model2[3]),
  std_error = c(standard_errors_model1[2], standard_errors_model1[3],
                standard_errors_model2[2], standard_errors_model2[3])
)

# Calculate 95% confidence intervals
coefs_ci_manual$conf_low <- coefs_ci_manual$estimate - 1.96 * coefs_ci_manual$std_error
coefs_ci_manual$conf_high <- coefs_ci_manual$estimate + 1.96 * coefs_ci_manual$std_error

# Set viridis color palette
viridis_palette <- viridis(2, option = "D", begin = 0.1, end = 0.7)

# Update term labels
coefs_ci_manual$term <- ifelse(coefs_ci_manual$term == "Geographic Distance", "Geographic Distance", coefs_ci_manual$term)
coefs_ci_manual$term <- ifelse(coefs_ci_manual$term == "Time Distance", "Time", coefs_ci_manual$term)

coefs_ci_manual$model <- factor(coefs_ci_manual$model, levels = c("Model 1", "Model 2"),
                                labels = c("Linear Model", "Linear Model, With Interaction"))

# Create a ggplot for coefficients and confidence intervals
coefficients_plot_DK <- ggplot(coefs_ci_manual, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = "Molecular Change", x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE)





# Urban vs countryside ----------
sampled_data_byzone<-readRDS(file="sampled_data_byzone.rds")
cophenetic_distances_byzone<-readRDS(file="cophenetic_distances_byzone.rds")
geographic_distance_matrix_byzone<-readRDS(file="geographic_distance_matrix_byzone.rds")
time_distance_matrix_byzone <- readRDS(file="time_distance_matrix_byzone.rds")

lower_tri_geographic_byzone <-  as.vector(lower_triangular(geographic_distance_matrix_byzone) / 10000)
lower_tri_time_byzone <- as.vector(abs(lower_triangular(time_distance_matrix_byzone)) / 7)
lower_tri_cophenetic_byzone <-  as.vector(lower_triangular(cophenetic_distances_byzone) * 29891)

combined_data_byzone <- data.frame(
  lower_tri_geographic_byzone = lower_tri_geographic_byzone,
  lower_tri_time_byzone = lower_tri_time_byzone,
  lower_tri_cophenetic_byzone = lower_tri_cophenetic_byzone)

geography_genetic_time_model_byzone <- biglm(lower_tri_cophenetic_byzone ~ lower_tri_geographic_byzone + lower_tri_time_byzone, data = combined_data_byzone)
summary(geography_genetic_time_model_byzone)





# Countryside / landzone
sampled_data_landzone<-readRDS(file="sampled_data_landzone.rds")
cophenetic_distances_landzone<-readRDS(file="cophenetic_distances_landzone.rds")
geographic_distance_matrix_landzone<-readRDS(file="geographic_distance_matrix_landzone.rds")
time_distance_matrix_landzone <- readRDS(file="time_distance_matrix_landzone.rds")

lower_tri_geographic_landzone <-  as.vector(lower_triangular(geographic_distance_matrix_landzone) / 10000)
lower_tri_time_landzone <- as.vector(abs(lower_triangular(time_distance_matrix_landzone)) / 7)
lower_tri_cophenetic_landzone <-  as.vector(lower_triangular(cophenetic_distances_landzone) * 29891)

combined_data_landzone <- data.frame(
  lower_tri_geographic_landzone = lower_tri_geographic_landzone,
  lower_tri_time_landzone = lower_tri_time_landzone,
  lower_tri_cophenetic_landzone = lower_tri_cophenetic_landzone)

geography_genetic_time_model_landzone <- biglm(lower_tri_cophenetic_landzone ~ lower_tri_geographic_landzone + lower_tri_time_landzone, data = combined_data_landzone)
summary(geography_genetic_time_model_landzone)















# Regression by region, 10,000 samples for each -------

# Nordjylland-------------
time_distance_matrix_Nordjylland <- readRDS(file="time_distance_matrix_Region_Nordjylland.rds")
geographic_distance_matrix_Nordjylland <- readRDS(file = "geographic_distance_matrix_Region_Nordjylland.rds")
cophenetic_distances_Nordjylland <- readRDS(file="cophenetic_distances_Region_Nordjylland.rds")
# Lower triangular + flattening
geographic_distance_matrix_parallel_vector_Nordjylland <- as.vector(lower_triangular(geographic_distance_matrix_Nordjylland) / 10000)
time_distance_matrix_Nordjylland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Nordjylland)) / 7)
cophenetic_distances_vector_Nordjylland <- as.vector(lower_triangular(cophenetic_distances_Nordjylland) * 29891)
combined_data_Nordjylland <- data.frame(
  time_distance_matrix_Nordjylland_vector = time_distance_matrix_Nordjylland_vector,
  cophenetic_distances_vector_Nordjylland = cophenetic_distances_vector_Nordjylland,
  geographic_distance_matrix_parallel_vector_Nordjylland = geographic_distance_matrix_parallel_vector_Nordjylland)
geography_genetic_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ geographic_distance_matrix_parallel_vector_Nordjylland + time_distance_matrix_Nordjylland_vector, data = combined_data_Nordjylland)
geography_genetic_model_interaction_model_Nordjylland <- biglm(cophenetic_distances_vector_Nordjylland ~ geographic_distance_matrix_parallel_vector_Nordjylland*time_distance_matrix_Nordjylland_vector, data = combined_data_Nordjylland)
#Plot:
# Get coefficients and their standard errors for the first model
coefficients_model1_Nordjylland <- coef(geography_genetic_model_interaction_model_Nordjylland)
standard_errors_model1_Nordjylland <- sqrt(diag(vcov(geography_genetic_model_interaction_model_Nordjylland)))
# Get coefficients and their standard errors for the second model
coefficients_model2_Nordjylland <- coef(geography_genetic_time_interaction_model)
standard_errors_model2_Nordjylland <- sqrt(diag(vcov(geography_genetic_time_interaction_model)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_Nordjylland <- data.frame(
  model = rep(c("Linear Model", "Linear Model, With Interaction"), each = 2),
  term = c("Geographic Distance", "Time"),
  estimate = c(coefficients_model1_Nordjylland[2], coefficients_model1_Nordjylland[3],
               coefficients_model2_Nordjylland[2], coefficients_model2_Nordjylland[3]),
  std_error = c(standard_errors_model1_Nordjylland[2], standard_errors_model1_Nordjylland[3],
                standard_errors_model2_Nordjylland[2], standard_errors_model2_Nordjylland[3])
)
coefs_ci_manual_Nordjylland$conf_low <- coefs_ci_manual_Nordjylland$estimate - 1.96 * coefs_ci_manual_Nordjylland$std_error
coefs_ci_manual_Nordjylland$conf_high <- coefs_ci_manual_Nordjylland$estimate + 1.96 * coefs_ci_manual_Nordjylland$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_Nordjylland <- ggplot(coefs_ci_manual_Nordjylland, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = NULL, x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE) +  # Remove legend for color and shape
  ggtitle(paste("Nordjylland"))
coefficients_plot_Nordjylland












# Midtjylland  -----
time_distance_matrix_Midtjylland <- readRDS(file="time_distance_matrix_Region_Midtjylland.rds")
geographic_distance_matrix_Midtjylland <- readRDS(file = "geographic_distance_matrix_Region_Midtjylland.rds")
cophenetic_distances_Midtjylland <- readRDS(file="cophenetic_distances_Region_Midtjylland.rds")
# Lower triangular + flattening
geographic_distance_matrix_parallel_vector_Midtjylland <- as.vector(lower_triangular(geographic_distance_matrix_Midtjylland) / 10000)
time_distance_matrix_Midtjylland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Midtjylland)) / 7)
cophenetic_distances_vector_Midtjylland <- as.vector(lower_triangular(cophenetic_distances_Midtjylland) * 29891)
combined_data_Midtjylland <- data.frame(
  time_distance_matrix_Midtjylland_vector = time_distance_matrix_Midtjylland_vector,
  cophenetic_distances_vector_Midtjylland = cophenetic_distances_vector_Midtjylland,
  geographic_distance_matrix_parallel_vector_Midtjylland = geographic_distance_matrix_parallel_vector_Midtjylland)
geography_genetic_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ geographic_distance_matrix_parallel_vector_Midtjylland + time_distance_matrix_Midtjylland_vector, data = combined_data_Midtjylland)
geography_genetic_model_interaction_model_Midtjylland <- biglm(cophenetic_distances_vector_Midtjylland ~ geographic_distance_matrix_parallel_vector_Midtjylland*time_distance_matrix_Midtjylland_vector, data = combined_data_Midtjylland)
#Plot:
# Get coefficients and their standard errors for the first model
coefficients_model1_Midtjylland <- coef(geography_genetic_model_Midtjylland)
standard_errors_model1_Midtjylland <- sqrt(diag(vcov(geography_genetic_model_Midtjylland)))
# Get coefficients and their standard errors for the second model
coefficients_model2_Midtjylland <- coef(geography_genetic_model_interaction_model_Midtjylland)
standard_errors_model2_Midtjylland <- sqrt(diag(vcov(geography_genetic_model_interaction_model_Midtjylland)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_Midtjylland <- data.frame(
  model = rep(c("Linear Model", "Linear Model, With Interaction"), each = 2),
  term = c("Geographic Distance", "Time"),
  estimate = c(coefficients_model1_Midtjylland[2], coefficients_model1_Midtjylland[3],
               coefficients_model2_Midtjylland[2], coefficients_model2_Midtjylland[3]),
  std_error = c(standard_errors_model1_Midtjylland[2], standard_errors_model1_Midtjylland[3],
                standard_errors_model2_Midtjylland[2], standard_errors_model2_Midtjylland[3])
)
coefs_ci_manual_Midtjylland$conf_low <- coefs_ci_manual_Midtjylland$estimate - 1.96 * coefs_ci_manual_Midtjylland$std_error
coefs_ci_manual_Midtjylland$conf_high <- coefs_ci_manual_Midtjylland$estimate + 1.96 * coefs_ci_manual_Midtjylland$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_Midtjylland <- ggplot(coefs_ci_manual_Midtjylland, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = NULL, x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE) +  # Remove legend for color and shape
  ggtitle(paste("Midtjylland"))
coefficients_plot_Midtjylland






# Syddanmark -------
time_distance_matrix_Syddanmark <- readRDS(file="time_distance_matrix_Region_Syddanmark.rds")
geographic_distance_matrix_Syddanmark <- readRDS(file = "geographic_distance_matrix_Region_Syddanmark.rds")
cophenetic_distances_Syddanmark <- readRDS(file="cophenetic_distances_Region_Syddanmark.rds")
# Lower triangular + flattening
geographic_distance_matrix_parallel_vector_Syddanmark <- as.vector(lower_triangular(geographic_distance_matrix_Syddanmark) / 10000)
time_distance_matrix_Syddanmark_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Syddanmark)) / 7)
cophenetic_distances_vector_Syddanmark <- as.vector(lower_triangular(cophenetic_distances_Syddanmark) * 29891)
combined_data_Syddanmark <- data.frame(
  time_distance_matrix_Syddanmark_vector = time_distance_matrix_Syddanmark_vector,
  cophenetic_distances_vector_Syddanmark = cophenetic_distances_vector_Syddanmark,
  geographic_distance_matrix_parallel_vector_Syddanmark = geographic_distance_matrix_parallel_vector_Syddanmark)
geography_genetic_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ geographic_distance_matrix_parallel_vector_Syddanmark + time_distance_matrix_Syddanmark_vector, data = combined_data_Syddanmark)
geography_genetic_model_interaction_model_Syddanmark <- biglm(cophenetic_distances_vector_Syddanmark ~ geographic_distance_matrix_parallel_vector_Syddanmark*time_distance_matrix_Syddanmark_vector, data = combined_data_Syddanmark)
#Plot:
# Get coefficients and their standard errors for the first model
coefficients_model1_Syddanmark <- coef(geography_genetic_model_Syddanmark)
standard_errors_model1_Syddanmark <- sqrt(diag(vcov(geography_genetic_model_Syddanmark)))
# Get coefficients and their standard errors for the second model
coefficients_model2_Syddanmark <- coef(geography_genetic_model_interaction_model_Syddanmark)
standard_errors_model2_Syddanmark <- sqrt(diag(vcov(geography_genetic_model_interaction_model_Syddanmark)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_Syddanmark <- data.frame(
  model = rep(c("Linear Model", "Linear Model, With Interaction"), each = 2),
  term = c("Geographic Distance", "Time"),
  estimate = c(coefficients_model1_Syddanmark[2], coefficients_model1_Syddanmark[3],
               coefficients_model2_Syddanmark[2], coefficients_model2_Syddanmark[3]),
  std_error = c(standard_errors_model1_Syddanmark[2], standard_errors_model1_Syddanmark[3],
                standard_errors_model2_Syddanmark[2], standard_errors_model2_Syddanmark[3])
)
coefs_ci_manual_Syddanmark$conf_low <- coefs_ci_manual_Syddanmark$estimate - 1.96 * coefs_ci_manual_Syddanmark$std_error
coefs_ci_manual_Syddanmark$conf_high <- coefs_ci_manual_Syddanmark$estimate + 1.96 * coefs_ci_manual_Syddanmark$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_Syddanmark <- ggplot(coefs_ci_manual_Syddanmark, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = NULL, x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE) +  # Remove legend for color and shape
  ggtitle(paste("Syddanmark"))
coefficients_plot_Syddanmark


# Hovedstaden -------
time_distance_matrix_Hovedstaden <- readRDS(file="time_distance_matrix_Region_Hovedstaden.rds")
geographic_distance_matrix_Hovedstaden <- readRDS(file = "geographic_distance_matrix_Region_Hovedstaden.rds")
cophenetic_distances_Hovedstaden <- readRDS(file="cophenetic_distances_Region_Hovedstaden.rds")
# Lower triangular + flattening
geographic_distance_matrix_parallel_vector_Hovedstaden <- as.vector(lower_triangular(geographic_distance_matrix_Hovedstaden) / 10000)
time_distance_matrix_Hovedstaden_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Hovedstaden)) / 7)
cophenetic_distances_vector_Hovedstaden <- as.vector(lower_triangular(cophenetic_distances_Hovedstaden) * 29891)
combined_data_Hovedstaden <- data.frame(
  time_distance_matrix_Hovedstaden_vector = time_distance_matrix_Hovedstaden_vector,
  cophenetic_distances_vector_Hovedstaden = cophenetic_distances_vector_Hovedstaden,
  geographic_distance_matrix_parallel_vector_Hovedstaden = geographic_distance_matrix_parallel_vector_Hovedstaden)
geography_genetic_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ geographic_distance_matrix_parallel_vector_Hovedstaden + time_distance_matrix_Hovedstaden_vector, data = combined_data_Hovedstaden)
geography_genetic_model_interaction_model_Hovedstaden <- biglm(cophenetic_distances_vector_Hovedstaden ~ geographic_distance_matrix_parallel_vector_Hovedstaden*time_distance_matrix_Hovedstaden_vector, data = combined_data_Hovedstaden)
#Plot:
# Get coefficients and their standard errors for the first model
coefficients_model1_Hovedstaden <- coef(geography_genetic_model_Hovedstaden)
standard_errors_model1_Hovedstaden <- sqrt(diag(vcov(geography_genetic_model_Hovedstaden)))
# Get coefficients and their standard errors for the second model
coefficients_model2_Hovedstaden <- coef(geography_genetic_model_interaction_model_Hovedstaden)
standard_errors_model2_Hovedstaden <- sqrt(diag(vcov(geography_genetic_model_interaction_model_Hovedstaden)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_Hovedstaden <- data.frame(
  model = rep(c("Linear Model", "Linear Model, With Interaction"), each = 2),
  term = c("Geographic Distance", "Time"),
  estimate = c(coefficients_model1_Hovedstaden[2], coefficients_model1_Hovedstaden[3],
               coefficients_model2_Hovedstaden[2], coefficients_model2_Hovedstaden[3]),
  std_error = c(standard_errors_model1_Hovedstaden[2], standard_errors_model1_Hovedstaden[3],
                standard_errors_model2_Hovedstaden[2], standard_errors_model2_Hovedstaden[3])
)
coefs_ci_manual_Hovedstaden$conf_low <- coefs_ci_manual_Hovedstaden$estimate - 1.96 * coefs_ci_manual_Hovedstaden$std_error
coefs_ci_manual_Hovedstaden$conf_high <- coefs_ci_manual_Hovedstaden$estimate + 1.96 * coefs_ci_manual_Hovedstaden$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_Hovedstaden <- ggplot(coefs_ci_manual_Hovedstaden, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = NULL, x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE) +  # Remove legend for color and shape
  ggtitle(paste("Hovedstaden"))
coefficients_plot_Hovedstaden


# Sjælland -------
time_distance_matrix_Sjælland <- readRDS(file="time_distance_matrix_Region_Sjælland.rds")
geographic_distance_matrix_Sjælland <- readRDS(file = "geographic_distance_matrix_Region_Sjælland.rds")
cophenetic_distances_Sjælland <- readRDS(file="cophenetic_distances_Region_Sjælland.rds")
# Lower triangular + flattening
geographic_distance_matrix_parallel_vector_Sjælland <- as.vector(lower_triangular(geographic_distance_matrix_Sjælland) / 10000)
time_distance_matrix_Sjælland_vector <- as.vector(abs(lower_triangular(time_distance_matrix_Sjælland)) / 7)
cophenetic_distances_vector_Sjælland <- as.vector(lower_triangular(cophenetic_distances_Sjælland) * 29891)
combined_data_Sjælland <- data.frame(
  time_distance_matrix_Sjælland_vector = time_distance_matrix_Sjælland_vector,
  cophenetic_distances_vector_Sjælland = cophenetic_distances_vector_Sjælland,
  geographic_distance_matrix_parallel_vector_Sjælland = geographic_distance_matrix_parallel_vector_Sjælland)
geography_genetic_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ geographic_distance_matrix_parallel_vector_Sjælland + time_distance_matrix_Sjælland_vector, data = combined_data_Sjælland)
geography_genetic_model_interaction_model_Sjælland <- biglm(cophenetic_distances_vector_Sjælland ~ geographic_distance_matrix_parallel_vector_Sjælland*time_distance_matrix_Sjælland_vector, data = combined_data_Sjælland)
#Plot:
# Get coefficients and their standard errors for the first model
coefficients_model1_Sjælland <- coef(geography_genetic_model_Sjælland)
standard_errors_model1_Sjælland <- sqrt(diag(vcov(geography_genetic_model_Sjælland)))
# Get coefficients and their standard errors for the second model
coefficients_model2_Sjælland <- coef(geography_genetic_model_interaction_model_Sjælland)
standard_errors_model2_Sjælland <- sqrt(diag(vcov(geography_genetic_model_interaction_model_Sjælland)))
# Combine coefficients and standard errors into a data frame
coefs_ci_manual_Sjælland <- data.frame(
  model = rep(c("Linear Model", "Linear Model, With Interaction"), each = 2),
  term = c("Geographic Distance", "Time"),
  estimate = c(coefficients_model1_Sjælland[2], coefficients_model1_Sjælland[3],
               coefficients_model2_Sjælland[2], coefficients_model2_Sjælland[3]),
  std_error = c(standard_errors_model1_Sjælland[2], standard_errors_model1_Sjælland[3],
                standard_errors_model2_Sjælland[2], standard_errors_model2_Sjælland[3])
)
coefs_ci_manual_Sjælland$conf_low <- coefs_ci_manual_Sjælland$estimate - 1.96 * coefs_ci_manual_Sjælland$std_error
coefs_ci_manual_Sjælland$conf_high <- coefs_ci_manual_Sjælland$estimate + 1.96 * coefs_ci_manual_Sjælland$std_error
# Create a ggplot for coefficients and confidence intervals
coefficients_plot_Sjælland <- ggplot(coefs_ci_manual_Sjælland, aes(x = model, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = NULL, x=NULL) +
  theme_minimal() +
  scale_color_manual(values = viridis_palette) +
  scale_shape_manual(values = c(16, 17)) +
  facet_wrap(~term, scales = "free_y", ncol = 1) +
  guides(color = FALSE, shape = FALSE) +  # Remove legend for color and shape
  ggtitle(paste("Sjælland"))
coefficients_plot_Sjælland




# New plot, interaction model, stratifying by region ----------
full_linear_interaction_model_results <- data.frame(
  Region = c("Denmark", "Denmark", "Nordjylland", "Nordjylland", "Syddanmark", "Syddanmark",
              "Hovedstaden", "Hovedstaden", "Midtjylland", "Midtjylland", "Sjælland", "Sjælland"),
  term = c("Geographic Distance", "Time", "Geographic Distance", "Time", "Geographic Distance", "Time",
           "Geographic Distance", "Time", "Geographic Distance", "Time", "Geographic Distance", "Time"),
  estimate = c(coefficients_model2[2], coefficients_model2[3],
               coefficients_model2_Nordjylland[2], coefficients_model2_Nordjylland[3],
               coefficients_model2_Syddanmark[2], coefficients_model2_Syddanmark[3],
               coefficients_model2_Hovedstaden[2], coefficients_model2_Hovedstaden[3],
               coefficients_model2_Midtjylland[2], coefficients_model2_Midtjylland[3],
               coefficients_model2_Sjælland[2], coefficients_model2_Sjælland[3]),
  std_error = c(standard_errors_model2[2], standard_errors_model2[3],
                standard_errors_model2_Nordjylland[2], standard_errors_model2_Nordjylland[3],
                standard_errors_model2_Syddanmark[2], standard_errors_model2_Syddanmark[3],
                standard_errors_model2_Hovedstaden[2], standard_errors_model2_Hovedstaden[3],
                standard_errors_model2_Midtjylland[2], standard_errors_model2_Midtjylland[3],
                standard_errors_model2_Sjælland[2], standard_errors_model2_Sjælland[3])
)

full_linear_interaction_model_results$conf_low <- full_linear_interaction_model_results$estimate - 1.96 * full_linear_interaction_model_results$std_error
full_linear_interaction_model_results$conf_high <- full_linear_interaction_model_results$estimate + 1.96 * full_linear_interaction_model_results$std_error

interaction_model_plot <- ggplot(full_linear_interaction_model_results, aes(x = Region, y = estimate, color = Region)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(
    aes(ymin = conf_low, ymax = conf_high),
    position = position_dodge(width = 0.8),
    width = 0.5
  ) +
  facet_wrap(~ term + ., scales = "free_y", ncol = 1, strip.position = "top") +
  theme_minimal() +
  labs(x = "Region", y = "Estimate", color = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_viridis_d() +
  ggtitle("With Interaction")
interaction_model_plot
ggsave(file="interaction_model_plot.pdf", interaction_model_plot, width = 4, height = 6)


# New plot, non-interaction model, stratifying by region ----------
full_linear_no_interaction_model_results <- data.frame(
  Region = c("Denmark", "Denmark", "Nordjylland", "Nordjylland", "Syddanmark", "Syddanmark",
             "Hovedstaden", "Hovedstaden", "Midtjylland", "Midtjylland", "Sjælland", "Sjælland"),
  term = c("Geographic Distance", "Time", "Geographic Distance", "Time", "Geographic Distance", "Time",
           "Geographic Distance", "Time", "Geographic Distance", "Time", "Geographic Distance", "Time"),
  estimate = c(coefficients_model1[2], coefficients_model1[3],
               coefficients_model1_Nordjylland[2], coefficients_model1_Nordjylland[3],
               coefficients_model1_Syddanmark[2], coefficients_model1_Syddanmark[3],
               coefficients_model1_Hovedstaden[2], coefficients_model1_Hovedstaden[3],
               coefficients_model1_Midtjylland[2], coefficients_model1_Midtjylland[3],
               coefficients_model1_Sjælland[2], coefficients_model1_Sjælland[3]),
  std_error = c(standard_errors_model1[2], standard_errors_model1[3],
                standard_errors_model1_Nordjylland[2], standard_errors_model1_Nordjylland[3],
                standard_errors_model1_Syddanmark[2], standard_errors_model1_Syddanmark[3],
                standard_errors_model1_Hovedstaden[2], standard_errors_model1_Hovedstaden[3],
                standard_errors_model1_Midtjylland[2], standard_errors_model1_Midtjylland[3],
                standard_errors_model1_Sjælland[2], standard_errors_model1_Sjælland[3])
)

full_linear_no_interaction_model_results$conf_low <- full_linear_no_interaction_model_results$estimate - 1.96 * full_linear_no_interaction_model_results$std_error
full_linear_no_interaction_model_results$conf_high <- full_linear_no_interaction_model_results$estimate + 1.96 * full_linear_no_interaction_model_results$std_error

non_interaction_model_plot <- ggplot(full_linear_no_interaction_model_results, aes(x = Region, y = estimate, color = Region)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(
    aes(ymin = conf_low, ymax = conf_high),
    position = position_dodge(width = 0.8),
    width = 0.5
  ) +
  facet_wrap(~ term + ., scales = "free_y", ncol = 1, strip.position = "top") +
  theme_minimal() +
  labs(x = "Region", y = "Estimate", color = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_viridis_d() +
  ggtitle("Without Interaction")
non_interaction_model_plot
ggsave(file="non_interaction_model_plot.pdf", non_interaction_model_plot, width = 4, height = 6)


# Combined plot --------
overall_plot <- non_interaction_model_plot + interaction_model_plot +
  plot_layout(ncol = 2)
ggsave(file="linear_regression_coefficient_plot.pdf", overall_plot, width = 8, height = 6)








# Correlation geo distance and phylogenetic distance using castor instead -----
install.packages(package_file_path, repos = NULL, type = "source")
library(castor)
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
                                                              correlation_method="pearson",
                                                              max_phylodistance = Inf,
                                                              Npermutations = 1000,
                                                              alternative = "two-sided",
                                                              radius = 6371)


