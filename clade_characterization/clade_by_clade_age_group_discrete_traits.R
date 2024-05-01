# Ancestral state reconstruction using age groups

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages  ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, mapdata, sf,
  rnaturalearth, rnaturalearthdata, ggplotify, phangorn, igraph, networkD3, heatmaply,
  pheatmap, RColorBrewer
)

cividis_palette <- viridisLite::cividis


# Load metadata -----
sequenced_individual_detailed_metadata <- readRDS(file = "")
region_name_mapping <- c(
  "1081" = "Nordjylland",
  "1082" = "Midtjylland",
  "1083" = "Syddanmark",
  "1084" = "Hovedstaden",
  "1085" = "Sjælland"
)
sequenced_individual_detailed_metadata$region <- factor(region_name_mapping[as.character(sequenced_individual_detailed_metadata$REGIONSKODE)])
all_positive_individuals <- readRDS(file = "")
IAR_denmark_data <- read.csv(file = "")
sequenced_individuals <- read.csv(file = "")

# Region colour mapping:
region_color_mapping <- data.frame(
  region = c("Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark"),
  color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

node_region_color_mapping <- data.frame(
  region = c("Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark"),
  color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

# Complete Phylogeny ---------
current_MAPLE_tree <- ape::read.tree("")
current_metadata <- read.csv(file = "")
selected_data <- sequenced_individual_detailed_metadata[, c("strain", "region", "age_at_infection", "REGIONSKODE")]
merged_data <- merge(current_metadata, selected_data, by = "strain", all.x = TRUE)
# Filter the tree to keep only tips present in merged_data
tips_to_keep <- sequenced_individual_detailed_metadata$strain
current_MAPLE_tree_filtered <- keep.tip(current_MAPLE_tree, tips_to_keep)

# Filter metadata to match IDs in the tree
filtered_data <- merged_data[merged_data$strain %in% current_MAPLE_tree_filtered$tip.label, ]

DTA_data <- filtered_data[, c("strain", "region", "age_at_infection", "REGIONSKODE")]
DTA_data_matched <- DTA_data[match(current_MAPLE_tree_filtered$tip.label, DTA_data$strain), ]
DTA_data_matched$age_group <- cut(DTA_data_matched$age_at_infection,
  breaks = c(0, 15, 30, 45, 60, 75, Inf),
  labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"),
  include.lowest = TRUE
)

states_age <- DTA_data_matched$age_group
states_region <- DTA_data_matched$region

# DTA
ace_output_age <- ace(states_age, phy = current_MAPLE_tree_filtered, type = "discrete", method = "ML")
ace_output_region <- ace(states_region, phy = current_MAPLE_tree_filtered, type = "discrete", method = "ML")

# Getting most likely states
i_age <- apply(ace_output_age$lik.anc, 1, which.max)

tip_age_trait <- data.frame(
  label = current_MAPLE_tree_filtered$tip.label,
  age_group = factor(DTA_data_matched$age_group[match(current_MAPLE_tree_filtered$tip.label, DTA_data_matched$strain)])
)

# Extract unique regions from the data
node_age_trait <- data.frame(
  label = current_MAPLE_tree_filtered$node.label,
  node_region = factor(i_age),
  age_group = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")[i_age]
)

# Create heatmap using ggplot2
# Creating transition matrix between all nodes
node_age_trait <- node_age_trait[, c("label", "age_group")]
combined_trait <- rbind(tip_age_trait, node_age_trait)
age_groups <- unique(combined_trait$age_group)
age_groups <- age_groups[!is.na(age_groups)]
transition_matrix_age <- matrix(0, nrow = length(age_groups), ncol = length(age_groups), dimnames = list(age_groups, age_groups))
correct_order <- c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")
transition_matrix_age <- transition_matrix_age[correct_order, correct_order]

for (i in 1:nrow(current_MAPLE_tree_filtered$edge)) {
  from_label <- as.numeric(rownames(combined_trait)[current_MAPLE_tree_filtered$edge[i, 1]])
  to_label <- as.numeric(rownames(combined_trait)[current_MAPLE_tree_filtered$edge[i, 2]])

  from_age <- combined_trait$age_group[from_label]
  to_age <- combined_trait$age_group[to_label]
  transition_matrix_age[from_age, to_age] <- transition_matrix_age[from_age, to_age] + 1
}

# Convert matrix to dataframe
transition_matrix_age <- as.matrix(transition_matrix_age)
valid_names <- c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")
rows_to_keep <- intersect(valid_names, rownames(transition_matrix_age))
cols_to_keep <- intersect(valid_names, colnames(transition_matrix_age))
transition_matrix_age <- transition_matrix_age[rows_to_keep, cols_to_keep, drop = FALSE]

df <- as.data.frame.table(transition_matrix_age)
df$Freq <- factor(df$Freq)
df$Freq <- as.numeric(as.character(df$Freq))
sorted_var1 <- sort(unique(df$Var1))
sorted_var2 <- sort(unique(df$Var2))
df$Var1 <- factor(df$Var1, levels = sorted_var1)
df$Var2 <- factor(df$Var2, levels = sorted_var2)

# Plot
full_tree_heatmap <- ggplot(df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c() + # Using continuous color scale
  theme_minimal() +
  labs(title = "Full Tree", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

full_tree_heatmap_text <- ggplot(df, aes(x = Var1, y = Var2, fill = Freq, label = Freq)) +
  geom_tile(color = "black") +
  geom_text(color = "black") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Full Tree", x = NULL, y = NULL) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

full_tree_heatmap_normalized <- as.ggplot(pheatmap(transition_matrix_age,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cividis(20),
  scale = "column",
  border_color = "black",
  fontsize = 8,
  display_numbers = TRUE,
  number_color = "black",
  legend = FALSE,
  main = "Full Tree",
  angle_col = 45
))


combined_plot <- ggarrange(full_tree_heatmap_text, full_tree_heatmap_normalized, ncol = 2)
ggsave("", combined_plot, width = 12, height = 6)

# Proportion test comparing diagonal vs. off-diagonal values --------
within_same_age <- diag(transition_matrix_age)
total_transmissions <- sum(transition_matrix_age)
total_within_same_age <- sum(within_same_age)
total_between_age <- total_transmissions - total_within_same_age
proportion_within_same_age <- total_within_same_age / total_transmissions
proportion_between_age <- total_between_age / total_transmissions
test_result <- prop.test(
  x = c(total_within_same_age, total_between_age),
  n = c(total_transmissions, total_transmissions),
  alternative = "two.sided",
  correct = FALSE
)
print(test_result)
