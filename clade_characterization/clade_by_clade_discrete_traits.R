# Script for clade-by-clade ancestral state reconstruction

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages ----------
library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, treeio, gridExtra, cowplot, mapdata, sf,
  rnaturalearth, rnaturalearthdata, ggplotify, phangorn, igraph, networkD3, heatmaply,
  pheatmap, RColorBrewer
)

cividis_palette <- viridisLite::cividis



# Load metadata -----
sequenced_individual_detailed_metadata <- readRDS(file = "")
region_name_mapping <- c(
  "1081" = "N",
  "1082" = "M",
  "1083" = "SY",
  "1084" = "H",
  "1085" = "SJ"
)
sequenced_individual_detailed_metadata$region <- factor(region_name_mapping[as.character(sequenced_individual_detailed_metadata$REGIONSKODE)])
all_positive_individuals <- readRDS(file = "")
IAR_denmark_data <- read.csv(file = "")
sequenced_individuals <- read.csv(file = "")

# Region colour mapping:
region_color_mapping <- data.frame(
  region = c("H", "M", "N", "SJ", "SY"),
  color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

node_region_color_mapping <- data.frame(
  region = c("H", "M", "N", "SJ", "SY"),
  color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

valid_names <- c("H", "M", "N", "SJ", "SY")




# Reading in trees for each clade ------

# Sample code for each clade; done for each clade --------
current_BEAST_tree <- ape::read.tree("")
current_metadata <- read.csv("")
merged_data <- merge(current_metadata, sequenced_individual_detailed_metadata, by = "PERSON_ID", all.x = TRUE, suffixes = c(".current", ".detailed"))
filtered_data <- merged_data[merged_data$strain.current %in% current_BEAST_tree$tip.label, ]
DTA_data <- filtered_data[, c("strain.current", "region")]
DTA_data_matched <- DTA_data[match(current_BEAST_tree$tip.label, DTA_data$strain.current), ]
states_region <- DTA_data_matched$region
# Set node labels
current_BEAST_tree$node.label <- 1:ape::Nnode(current_BEAST_tree)
# Run ace for nodes
ace_output_nodes <- ape::ace(states_region, phy = current_BEAST_tree, type = "discrete", method = "ML")
i_region_nodes <- apply(ace_output_nodes$lik.anc, 1, which.max)

tip_region_trait <- data.frame(
  label = current_BEAST_tree$tip.label,
  region = factor(filtered_data$region[match(current_BEAST_tree$tip.label, filtered_data$strain.current)])
)

# Extract unique regions from the data
node_region_trait <- data.frame(
  label = current_BEAST_tree$node.label,
  node_region = factor(i_region_nodes),
  region_name = c("H", "M", "N", "SJ", "SY")[i_region_nodes]
)

# Creating tree
mrsd <- as.Date(max(current_metadata$DateSamplingLinelist))
decimal_dates <- lubridate::decimal_date(as.Date(current_metadata$DateSamplingLinelist))
viridis_palette <- viridis(nlevels(factor(tip_region_trait$region)))
tree_gg <- ggtree(current_BEAST_tree, mrsd = mrsd) + theme_tree2()

num_tips <- length(current_BEAST_tree$tip.label)
subtitle <- paste("n =", num_tips)
tree_gg <- tree_gg +
  labs(subtitle = subtitle)

tree1 <- tree_gg %<+% tip_region_trait +
  geom_tippoint(aes(color = region), alpha = 0.75) +
  scale_color_manual(values = setNames(region_color_mapping$color, region_color_mapping$region)) +
  theme(legend.position = "none") +
  ggtitle("")

tree <- tree1 %<+% node_region_trait +
  geom_nodepoint(aes(color = region_name, alpha = 0.75)) +
  scale_color_manual(values = setNames(node_region_color_mapping$color, node_region_color_mapping$region))

# Create heatmap using ggplot2
# Creating transition matrix between all nodes
names(tip_region_trait)[names(tip_region_trait) == "region"] <- "region_name"
node_region_trait <- node_region_trait[, c("label", "region_name")]
combined_trait <- rbind(tip_region_trait, node_region_trait)
transition_matrix_region <- matrix(0, nrow = length(valid_names), ncol = length(valid_names), dimnames = list(valid_names, valid_names))
for (i in 1:nrow(current_BEAST_tree$edge)) {
  from_label <- as.numeric(rownames(combined_trait)[current_BEAST_tree$edge[i, 1]])
  to_label <- as.numeric(rownames(combined_trait)[current_BEAST_tree$edge[i, 2]])

  from_region <- combined_trait$region_name[from_label]
  to_region <- combined_trait$region_name[to_label]

  # Check for NA values
  if (!is.na(from_region) && !is.na(to_region)) {
    transition_matrix_region[as.character(from_region), as.character(to_region)] <- transition_matrix_region[as.character(from_region), as.character(to_region)] + 1
  }
}


# Convert matrix to dataframe
transition_matrix_region <- as.matrix(transition_matrix_region)

df <- as.data.frame.table(transition_matrix_region)
df$Freq <- factor(df$Freq)
df$Freq <- as.numeric(as.character(df$Freq))
sorted_var1 <- sort(unique(df$Var1))
sorted_var2 <- sort(unique(df$Var2))
df$Var1 <- factor(df$Var1, levels = sorted_var1)
df$Var2 <- factor(df$Var2, levels = sorted_var2)

# Plot the heatmap
tree_heatmap <- ggplot(df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tree_heatmap_text <- ggplot(df, aes(x = Var1, y = Var2, fill = Freq, label = Freq)) +
  geom_tile(color = "black") +
  geom_text(color = "black") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "", x = NULL, y = NULL) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tree_heatmap_normalized <- as.ggplot(pheatmap(transition_matrix_region,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cividis(20),
  scale = "column",
  border_color = "black",
  fontsize = 12,
  display_numbers = FALSE,
  number_color = "black",
  legend = FALSE,
  main = "",
  angle_col = 45
))





# Complete Phylogeny/Tree ---------
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
i_region <- apply(ace_output_region$lik.anc, 1, which.max)
i_age <- apply(ace_output_age$lik.anc, 1, which.max)

tip_region_trait <- data.frame(
  label = current_MAPLE_tree_filtered$tip.label,
  region = factor(DTA_data_matched$region[match(current_MAPLE_tree_filtered$tip.label, DTA_data_matched$strain)])
)

# Extract unique regions from the data
node_region_trait <- data.frame(
  label = current_MAPLE_tree_filtered$node.label,
  node_region = factor(i_region),
  region_name = c("H", "M", "N", "SJ", "SY")[i_region]
)

# Create heatmap using ggplot2
# Creating transition matrix between all nodes
names(tip_region_trait)[names(tip_region_trait) == "region"] <- "region_name"
node_region_trait <- node_region_trait[, c("label", "region_name")]
combined_trait <- rbind(tip_region_trait, node_region_trait)
regions <- unique(combined_trait$region_name)
regions <- regions[!is.na(regions)]
transition_matrix_region <- matrix(0, nrow = length(regions), ncol = length(regions), dimnames = list(regions, regions))
for (i in 1:nrow(current_MAPLE_tree$edge)) {
  from_label <- as.numeric(rownames(combined_trait)[current_MAPLE_tree$edge[i, 1]])
  to_label <- as.numeric(rownames(combined_trait)[current_MAPLE_tree$edge[i, 2]])

  from_region <- combined_trait$region_name[from_label]
  to_region <- combined_trait$region_name[to_label]

  # Check for NA values
  if (!is.na(from_region) && !is.na(to_region)) {
    transition_matrix_region[as.character(from_region), as.character(to_region)] <- transition_matrix_region[as.character(from_region), as.character(to_region)] + 1
  }
}


# Convert matrix to dataframe
transition_matrix_region <- as.matrix(transition_matrix_region)

df <- as.data.frame.table(transition_matrix_region)
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

full_tree_heatmap_normalized <- as.ggplot(pheatmap(transition_matrix_region,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cividis(20),
  scale = "column", # You can adjust scaling as needed
  border_color = "black",
  fontsize = 12,
  display_numbers = FALSE,
  number_color = "black",
  legend = FALSE,
  main = "Full Tree",
  angle_col = 45
))

# Getting legend
legend_plot <- as.ggplot(pheatmap(transition_matrix_region,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cividis(20),
  scale = "column",
  border_color = "black",
  fontsize = 12,
  display_numbers = FALSE,
  number_color = "black",
  legend = TRUE,
  angle_col = 45
))





# Extract legend ------
tree_gg <- ggtree(current_BEAST_tree, mrsd = mrsd) + theme_tree2()
legend_tree <- tree_gg %<+% tip_region_trait +
  geom_tippoint(aes(color = region), alpha = 0.75) +
  scale_color_manual(values = setNames(region_color_mapping$color, region_color_mapping$region)) +
  theme(legend.position = "none") +
  theme(legend.position = "bottom") +
  ggtitle("Omicron BA.1")

legend <- cowplot::get_legend(legend_tree)


# Combine all the trees with the legend --------
combined_DTA_trees_with_legend <- grid.arrange(, legend, ncol = 4)
ggsave("", combined_DTA_trees_with_legend, width = 15, height = 20)


# Heatmap Grid, z-scored --------
heatmap_grid_z_Score <- grid.arrange(, NULL,
  full_tree_heatmap_normalized,
  ncol = 3
)
ggsave("", heatmap_grid_z_Score, width = 15, height = 20)

# Heatmap Grid, colours and numbers --------
heatmap_grid <- grid.arrange(, full_tree_heatmap,
  legend,
  ncol = 4
)
ggsave("", tree_and_heatmap_grid, width = 15, height = 20)


# Heatmap with Numbers and Trees Grid --------
tree_numbers_and_heatmap_grid <- grid.arrange(full_tree_heatmap_text,
  legend,
  ncol = 4
)
ggsave("", tree_numbers_and_heatmap_grid, width = 15, height = 20)


# Heatmap with Numbers and Trees Grid --------
tree_z_score_heatmap_grid <- grid.arrange(, legend, ncol = 4)
ggsave("", tree_z_score_heatmap_grid, width = 15, height = 20)


# Final Heatmap with Fewer Trees --------
final_tree_z_score_heatmap_grid <- grid.arrange(,
  full_tree_heatmap_normalized,
  legend_plot,
  ncol = 4
)
ggsave("", final_tree_z_score_heatmap_grid, width = 15, height = 20)
