# Average cophenetic distance between people over time

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, ggplot, ggplotify,
  phangorn, heatmaply, RColorBrewer, graphics, purrr, future.apply, geosphere, patchwork, coefplot
)

# Data Preparation (Tree and Metadata) ------------------------------
sequenced_individuals <- readRDS(file = "")
all_individual <- read.csv(file = "")
distance_tree <- read.tree("")

# Total sample -------------
unique_dates <- unique(sequenced_individuals$date)
dates_2021 <- unique_dates[
  as.Date(unique_dates) >= as.Date("2021-01-01") & as.Date(unique_dates) <= as.Date("2021-12-31")
]
dates_2021 <- sort(dates_2021)

# Initialize vectors to store results
mean_distances <- numeric(length(dates_2021))
median_distances <- numeric(length(dates_2021))
ci_lows <- numeric(length(dates_2021))
ci_highs <- numeric(length(dates_2021))

# Create a parallel backend
cores <- 8
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel loop through each date
foreach(i = seq_along(dates_2021), .combine = "c") %dopar% {
  date <- dates_2021[i]
  # Subset individuals by date
  individuals_on_date <- sequenced_individuals[strftime(sequenced_individuals$date, "%Y-%m-%d") == date, ]

  library(ape)

  # Subset tree by strains
  matched_strains <- match(distance_tree$tip.label, individuals_on_date$strain)
  tips_to_keep <- distance_tree$tip.label[!is.na(matched_strains)]
  pruned_tree <- ape::keep.tip(distance_tree, tips_to_keep)

  # Calculate cophenetic distances for relevant individuals
  relevant_cophenetic_distances <- cophenetic(pruned_tree)

  # Extract lower diagonal of the distance matrix
  lower_triangular <- relevant_cophenetic_distances[lower.tri(relevant_cophenetic_distances)]

  # Calculate mean, median, and confidence intervals
  mean_distance <- mean(lower_triangular)
  median_distance <- median(lower_triangular)
  ci <- t.test(lower_triangular)$conf.int
  ci_low <- ci[1]
  ci_high <- ci[2]

  # Return results
  c(mean_distance, median_distance, ci_low, ci_high)
} -> results

stopCluster(cl)

# Assign results to vectors
mean_distances <- results[seq(1, length(results), by = 4)]
median_distances <- results[seq(2, length(results), by = 4)]
ci_lows <- results[seq(3, length(results), by = 4)]
ci_highs <- results[seq(4, length(results), by = 4)]

saveRDS(mean_distances, file = "")
saveRDS(median_distances, file = "")
saveRDS(ci_lows, file = "")
saveRDS(ci_highs, file = "")

# Create a data frame for results
total_sample_data <- data.frame(
  date = as.Date(dates_2021),
  mean_distance = mean_distances,
  median_distance = median_distances,
  ci_low = ci_lows,
  ci_high = ci_highs
)


# Plot mean distances with 95% confidence bands
mean_plot <- ggplot(total_sample_data, aes(x = date, y = mean_distance)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "black", alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "red", linetype = "dashed") +
  labs(title = "Mean Cophenetic Distance over Time", x = "Date", y = "Mean Distance") +
  theme_minimal()



# Same process but for each region ----------
region_name_mapping <- c(
  "1081" = "Region Nordjylland",
  "1082" = "Region Midtjylland",
  "1083" = "Region Syddanmark",
  "1084" = "Region Hovedstaden",
  "1085" = "Region Sjælland"
)


# Function to perform analysis for each region
analyze_region <- function(region_data) {
  # Initialize vectors to store results
  mean_distances <- numeric(length(dates_2021))
  median_distances <- numeric(length(dates_2021))
  ci_lows <- numeric(length(dates_2021))
  ci_highs <- numeric(length(dates_2021))

  # Loop through each date
  for (i in seq_along(dates_2021)) {
    date <- dates_2021[i]
    # Subset individuals by date
    individuals_on_date <- region_data[strftime(region_data$date, "%Y-%m-%d") == date, ]

    library(ape)
    # Subset tree by strains
    matched_strains <- match(distance_tree$tip.label, individuals_on_date$strain)
    tips_to_keep <- distance_tree$tip.label[!is.na(matched_strains)]
    pruned_tree <- ape::keep.tip(distance_tree, tips_to_keep)

    # Calculate cophenetic distances for relevant individuals
    relevant_cophenetic_distances <- cophenetic(pruned_tree)

    # Extract lower diagonal of the distance matrix
    lower_triangular <- relevant_cophenetic_distances[lower.tri(relevant_cophenetic_distances)]

    # Calculate mean, median, and confidence intervals
    mean_distances[i] <- mean(lower_triangular)
    median_distances[i] <- median(lower_triangular)
    ci <- t.test(lower_triangular)$conf.int
    ci_lows[i] <- ci[1]
    ci_highs[i] <- ci[2]
  }

  # Create a data frame with results for the region
  region_results <- data.frame(
    date = as.Date(dates_2021),
    mean_distance = mean_distances,
    median_distance = median_distances,
    ci_low = ci_lows,
    ci_high = ci_highs
  )

  return(region_results)
}

# Register parallel backend
cl <- makeCluster(8)
registerDoParallel(cl)

filtered_data <- sequenced_individuals %>% filter(!is.na(REGIONSKODE))

# Split data by 'REGIONSKODE' and perform analysis for each region in parallel
results_by_region <- foreach(
  region_data = split(filtered_data, filtered_data$REGIONSKODE),
  .combine = rbind
) %dopar% {
  analyze_region(region_data)
}
# Stop the parallel backend
stopCluster(cl)

# Saving
results_by_region$Region <- rep(c(1081, 1082, 1083, 1084, 1085), each = 365)
results_by_region$region_name <- factor(region_name_mapping[as.character(results_by_region$Region)])
saveRDS(results_by_region, file = "results_by_region.rds")



# Merging -
total_sample_data$region_name <- "Full"
total_sample_data$Region <- 1
merged_data <- rbind(results_by_region, total_sample_data)

# Plotting ----
mean_plot <- ggplot(merged_data, aes(x = date, y = mean_distance)) +
  geom_line(aes(color = region_name)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = region_name), alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = "Mean Cophenetic Distance over Time",
    x = "Date",
    y = "Mean Distance",
    color = "Region"
  ) +
  theme_minimal() +
  facet_wrap(~region_name, scales = "free_y", ncol = 1)

ggplot(merged_data, aes(x = date, y = mean_distance * 29891, color = region_name)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "black", alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette
  labs(
    title = "Mean Cophenetic Distance over Time",
    x = "Date",
    y = "Mean Distance",
    color = "Region"
  ) +
  theme_minimal()

ggplot(merged_data, aes(x = date, y = median_distance * 29891, color = region_name)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "black", alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette
  labs(
    title = "Median Cophenetic Distance over Time",
    x = "Date",
    y = "Median Distance",
    color = "Region"
  ) +
  theme_minimal()


# Create a list of combinations of regions
region_combinations <- combn(unique(all_data$region_name), 2, simplify = TRUE)

# Function to create plots for each region combination
create_comparison_plot <- function(region1, region2) {
  # Filter data for the two regions
  comparison_data <- all_data %>%
    filter(region_name %in% c(region1, region2))

  # Plot mean distances for the two regions with viridis colors
  ggplot(comparison_data, aes(x = date, y = mean_distance, color = region_name)) +
    geom_line() +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "black", alpha = 0.3) +
    scale_color_viridis(discrete = TRUE) + # Use viridis color palette
    labs(title = paste("Mean Cophenetic Distance for", region1, "and", region2), x = "Date", y = "Mean Distance", color = "Region") +
    theme_minimal()
}

# Create plots for each region combination
comparison_plots <- lapply(seq_along(ncol(region_combinations)), function(i) {
  create_comparison_plot(region_combinations[1, i], region_combinations[2, i])
})

# Arrange plots in a grid
full_plot <- grid.arrange(grobs = comparison_plots, ncol = 2)
