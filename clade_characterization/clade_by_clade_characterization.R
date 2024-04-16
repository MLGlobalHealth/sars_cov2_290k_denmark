# Script for clade-by-clade description

library(
  ape, phytools, TreeTools, dplyr, tidyverse, data.table, dbplyr, lubridate,
  rlang, foreach, doParallel, DSTora, ROracle, DSTcolectica, DSTdb, DBI,
  parallel, ggsignif, viridis, ggtree, ggpubr, treeio, gridExtra, cowplot, mapdata, sf,
  rnaturalearth, rnaturalearthdata, ggplotify, phangorn, igraph, networkD3, heatmaply,
  pheatmap, RColorBrewer
)


# Establish connection
drv <- dbDriver("Oracle")
conn <- DSTora::OraGenvej("", dbuser = "")
conn2 <- DSTora::OraGenvej("", dbuser = "")
conn3 <- DSTora::OraGenvej("", dbuser = "")

# To access data tables

COVID_TEST <- dbReadTable(
  conn = conn,
  name = "",
  schema = ""
)

COVID_VAC <- dbReadTable(
  conn = conn,
  name = "",
  schema = ""
)


dbListTables(conn2, schema = "")


con <- get_conn(dbname = "")
dbListTables(con, schema = "")

lifelines <- dbReadTable(
  conn = conn,
  name = "",
  schema = ""
)
lifelines_koen <- dbReadTable(
  conn = conn,
  name = "",
  schema = ""
)

query <- ""
geo_address <- dbGetQuery(conn = con, statement = query)

general_address <- dbReadTable(
  conn = conn2,
  name = "",
  schema = ""
)



# All sequenced individuals:
sequenced_individuals <- read.csv(file = "")
all_positive_individuals <- COVID_TEST %>%
  filter(
    !is.na(PERSON_ID),
    between(PRDATE, as.Date("2021-01-01"), as.Date("2021-12-31"))
  ) %>%
  filter(SVARRESULTAT == 1)
all_positive_PCR_individuals <- all_positive_individuals %>%
  filter(CASEDEF == "SARS2") %>%
  filter(SVARRESULTAT == 1)

# Creating dataframe with metadata for sequenced individuals --------
common_person_ids <- sequenced_individuals %>%
  select(PERSON_ID)
# Step 1: Filter 'lifelines' for rows where 'PERSON_ID' is in 'sequenced_individuals'
filtered_lifelines <- lifelines %>%
  semi_join(common_person_ids, by = "PERSON_ID")
# Step 2: Join 'filtered_lifelines' with 'lifelines' on 'PERSON_ID' and take the first 'BIRTHDAY'
merged_data <- sequenced_individuals %>%
  left_join(
    filtered_lifelines %>%
      group_by(PERSON_ID) %>%
      summarize(BIRTHDAY = first(BIRTHDAY)),
    by = "PERSON_ID"
  )
# Step 3: Join 'merged_data' with 'lifelines_koen' on 'PERSON_ID'
sequenced_individual_detailed_metadata <- merged_data %>%
  left_join(
    lifelines_koen %>%
      group_by(PERSON_ID) %>%
      summarize(KOEN = first(KOEN)),
    by = "PERSON_ID"
  )
# Step 4: Join 'sequenced_individual_detailed_metadata' with 'geo_address' on 'PERSON_ID'
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  left_join(geo_address, by = "PERSON_ID")
# Step 5: Filter rows based on date range
# Assuming 'BOP_VFRA' and 'BOP_VTIL' are in POSIXct format
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  mutate(
    BOP_VFRA = as.Date(BOP_VFRA),
    BOP_VTIL = as.Date(BOP_VTIL)
  )
# Now, you can proceed with the filter
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  filter(date >= BOP_VFRA & date <= BOP_VTIL)
# Only keeping people that made the final tree
full_timetree <- read.tree("final_pruned_no_matOptimize_chronumental_consensus_2021_tree.tree")
tip_labels_tree <- full_timetree$tip.label
# Filter rows based on whether 'strain' is in tip labels
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  filter(strain %in% tip_labels_tree)
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  left_join(general_address, by = "ID")
# Creating new unique ADDRESS ID
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  mutate(UNIQUE_ADDRESS_ID = paste(KOM.x, VEJKODE.x, HUSNR.x, ETAGE, SIDEDOER, sep = "_"))
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  filter(!is.na(ID))
sequenced_individual_detailed_metadata$age_at_infection <-
  as.numeric(difftime(sequenced_individual_detailed_metadata$date,
    sequenced_individual_detailed_metadata$BIRTHDAY,
    units = "days"
  )) / 365.25
sequenced_individual_detailed_metadata$date <- as.Date(sequenced_individual_detailed_metadata$date, format = "%Y-%m-%d")



# Vaccination data:
COVID_VAC$FIRST_VACCINEDATE <- as.Date(COVID_VAC$FIRST_VACCINEDATE, format = "%Y-%m-%d")
COVID_VAC$SECOND_VACCINEDATE <- as.Date(COVID_VAC$SECOND_VACCINEDATE, format = "%Y-%m-%d")
sequenced_individual_detailed_metadata$date <- as.Date(sequenced_individual_detailed_metadata$date, format = "%Y-%m-%d")
sequenced_individual_detailed_metadata <- merge(sequenced_individual_detailed_metadata, COVID_VAC, by = "PERSON_ID", all.x = TRUE, suffixes = c(".seq", ".vac"))
sequenced_individual_detailed_metadata$Partially_vaccinated_at_exposure <- as.integer(ifelse(is.na(sequenced_individual_detailed_metadata$FIRST_VACCINEDATE), 0, sequenced_individual_detailed_metadata$FIRST_VACCINEDATE < sequenced_individual_detailed_metadata$date))
sequenced_individual_detailed_metadata$Fully_vaccinated_at_exposure <- as.integer(ifelse(is.na(sequenced_individual_detailed_metadata$SECOND_VACCINEDATE), 0, sequenced_individual_detailed_metadata$SECOND_VACCINEDATE < sequenced_individual_detailed_metadata$date))


# Load metadata -----
sequenced_individual_detailed_metadata <- readRDS(file = "")
region_name_mapping <- c(
  "1081" = "Region Nordjylland",
  "1082" = "Region Midtjylland",
  "1083" = "Region Syddanmark",
  "1084" = "Region Hovedstaden",
  "1085" = "Region Sjælland"
)
sequenced_individual_detailed_metadata$region <- factor(region_name_mapping[as.character(sequenced_individual_detailed_metadata$REGIONSKODE)])


all_positive_individuals <- readRDS(file = "")
IAR_denmark_data <- read.csv(file = "")
sequenced_individuals <- read.csv(file = "")

# Getting map of Denmark ------
world <- ne_countries(scale = "medium", returnclass = "sf")
# Remove Faroe Islands from the data
mainland_denmark <- subset(world, name_long == "Denmark" & iso_a2 != "FO")
denmark_bbox <- st_bbox(mainland_denmark)
# Create a ggplot with the specified bounding box
gg_denmark <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "white", size = 0.1) +
  geom_sf(data = mainland_denmark, fill = "white", color = "black", size = 0.5) +
  coord_sf(
    xlim = c(denmark_bbox$xmin, denmark_bbox$xmax),
    ylim = c(denmark_bbox$ymin, denmark_bbox$ymax)
  ) +
  theme_void()
print(gg_denmark)


# Functions ---------
# Check if any negative branch lengths
has_negative_lengths <- function(tree) {
  for (edge in tree$edge.length) {
    if (edge < 0) {
      return(TRUE)
    }
  }
  return(FALSE)
}
# Collapsing all negative branch lengths to zero - negative branch lengths is a known BEAST problem.
collapse_negative_lengths <- function(tree) {
  for (i in 1:length(tree$edge.length)) {
    if (tree$edge.length[i] < 0) {
      tree$edge.length[i] <- 0.00000001
    }
  }
  return(tree)
}

# Function to remove text after the last underscore
remove_after_last_underscore <- function(label) {
  gsub("_[^_]+$", "", label)
}

create_variant_map <- function(variant_name, tree, data, xlim, ylim, region_color_mapping) {
  # Create a data frame for tip metadata
  tip_metadata <- data.frame(
    label = tree$tip.label,
    lat = NA,
    long = NA,
    region = NA
  )

  # Region colour mapping:
  region_color_mapping <- data.frame(
    region = c("Region Hovedstaden", "Region Midtjylland", "Region Nordjylland", "Region Sjælland", "Region Syddanmark"),
    color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
  )


  # Match labels and indices
  matching_labels <- intersect(tree$tip.label, data$strain)
  matching_indices <- match(matching_labels, tree$tip.label)

  # Fill in metadata for matching labels
  tip_metadata[matching_indices, c("lat", "long", "region")] <- data[data$strain %in% matching_labels, c("latitude", "longitude", "region")]

  # Set row names
  rownames(tip_metadata) <- tip_metadata$label

  # Remove the 'label' column
  tip_metadata$label <- NULL

  # Generate phylo.to.map object
  danish_obj <- phylo.to.map(tree, tip_metadata, database = "worldHires", regions = "Denmark", xlim = xlim, ylim = ylim, plot = FALSE)

  # Use the specified color mapping for regions
  cols <- setNames(region_color_mapping$color[match(tip_metadata$region, region_color_mapping$region)], tree$tip.label)

  # Plot the map
  variant_tree_on_map <- plot(danish_obj, colors = cols, direction = "rightwards", ftype = "off", cex.points = c(0, 0.5), psize = 0.5, pts = FALSE, lwd = 0.1)

  return(variant_tree_on_map)
}

create_point_map <- function(data, xlim, ylim, region_color_mapping) {
  # Getting map of Denmark
  world <- ne_countries(scale = "medium", returnclass = "sf")
  mainland_denmark <- subset(world, name_long == "Denmark" & iso_a2 != "FO")
  denmark_bbox <- st_bbox(mainland_denmark)

  # Region colour mapping:
  region_color_mapping <- data.frame(
    region = c("Region Hovedstaden", "Region Midtjylland", "Region Nordjylland", "Region Sjælland", "Region Syddanmark"),
    color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
  )


  # Create a ggplot object with the specified bounding box
  gg_map <- ggplot() +
    geom_sf(data = world, fill = "lightgrey", color = "white", size = 0.1) +
    geom_sf(data = mainland_denmark, fill = "white", color = "black", size = 0.5) +
    coord_sf(
      xlim = c(denmark_bbox$xmin, denmark_bbox$xmax),
      ylim = c(denmark_bbox$ymin, denmark_bbox$ymax)
    ) +
    theme_void() +
    # Overlay points on the map
    geom_point(data = data, aes(x = longitude, y = latitude, color = region), size = 1.5) +
    scale_color_manual(values = setNames(region_color_mapping$color, region_color_mapping$region)) +
    xlim(xlim) +
    ylim(ylim) + # Set the plot limits
    theme_void() + # Set the plot limits
    guides(color = "none") # Remove the legend

  return(gg_map)
}





# Region colour mapping ---------
region_color_mapping <- data.frame(
  region = c("Region Hovedstaden", "Region Midtjylland", "Region Nordjylland", "Region Sjælland", "Region Syddanmark"),
  color = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)





# Reading in trees for each clade ------

# Code for each clade, repeated for each clade --------
current_BEAST_tree <- ape::read.nexus()
current_BEAST_tree <- collapse_negative_lengths(current_BEAST_tree)
current_BEAST_tree$tip.label <- sapply(current_BEAST_tree$tip.label, remove_after_last_underscore)
current_MAPLE_tree <- ape::read.tree("")
current_metadata <- read.csv("")
print(min(current_metadata$DateSamplingLinelist))
print(max(current_metadata$DateSamplingLinelist))
length(current_BEAST_tree$tip.label)
# Getting expected sampling fraction in this period:
filtered_IAR_data <- subset(IAR_denmark_data, date >= min(current_metadata$DateSamplingLinelist) & date <= max(current_metadata$DateSamplingLinelist))
average_iar <- mean(filtered_IAR_data$iar)
filtered_positive_individuals <- subset(all_positive_individuals, PRDATE >= min(current_metadata$DateSamplingLinelist) & PRDATE <= max(current_metadata$DateSamplingLinelist))
filtered_total_sequenced_individuals <- subset(sequenced_individuals, date >= min(current_metadata$DateSamplingLinelist) & date <= max(current_metadata$DateSamplingLinelist))
percentage_sequenced <- length(unique(filtered_total_sequenced_individuals$PERSON_ID)) / length(unique(filtered_positive_individuals$PERSON_ID))
expected_sampling_fraction <- average_iar * percentage_sequenced
print(expected_sampling_fraction) # 0.6210456
# Age
filtered_data <- sequenced_individual_detailed_metadata[sequenced_individual_detailed_metadata$PERSON_ID %in% current_metadata$PERSON_ID, ]
filtered_data <- filtered_data[filtered_data$strain %in% current_MAPLE_tree$tip.label, ]
summary(filtered_data$age_at_infection)
# Vaccination
table(filtered_data$Fully_vaccinated_at_exposure)
table(filtered_data$Partially_vaccinated_at_exposure)
# Region
table(filtered_data$region)
percentage_in_each_group <- prop.table(table(filtered_data$region)) * 100
print(percentage_in_each_group)
# Scorpio_call
table(filtered_data$scorpio_call)
percentage_in_each_group_scorpio <- prop.table(table(filtered_data$scorpio_call)) * 100
print(percentage_in_each_group_scorpio)
# Visualizing tree by region
tip_region_trait <- data.frame(
  label = current_BEAST_tree$tip.label,
  region = factor(sequenced_individual_detailed_metadata$region[match(current_BEAST_tree$tip.label, sequenced_individual_detailed_metadata$strain)])
)
mrsd <- as.Date(max(current_metadata$DateSamplingLinelist))
decimal_dates <- lubridate::decimal_date(as.Date(current_metadata$DateSamplingLinelist))
viridis_palette <- viridis(nlevels(factor(tip_region_trait$region)))
tree_gg <- ggtree(current_BEAST_tree, mrsd = mrsd) + theme_tree2()
Alpha_B.1.1_tree <- tree_gg %<+% tip_region_trait +
  geom_tippoint(aes(color = region), alpha = 0.75) +
  scale_color_manual(values = setNames(region_color_mapping$color, region_color_mapping$region)) +
  theme(legend.position = "none") +
  ggtitle("")
# Tree and map
tree_on_map_points <- create_point_map(filtered_data, xlim = c(7, 16), ylim = c(53, 60))
# Diffusion analysis
tip_labels_order <- current_BEAST_tree$tip.label
filtered_data <- filtered_data %>%
  filter(!is.na(latitude) & !is.na(longitude) & strain %in% tip_labels_order)
pruned_tree <- drop.tip(current_BEAST_tree, setdiff(tip_labels_order, filtered_data$strain))
new_tip_order <- pruned_tree$tip.label
tip_latitudes <- filtered_data$latitude[match(new_tip_order, filtered_data$strain)]
tip_longitudes <- filtered_data$longitude[match(new_tip_order, filtered_data$strain)]





# Combined Tree Figure ---------

# Extract legend from one of the figures
legend_plot <- tree_gg %<+% tip_region_trait +
  geom_tippoint(aes(color = region), alpha = 0.75) +
  scale_color_manual(values = viridis_palette) +
  labs(color = "Region")
legend_plot
legend <- cowplot::get_legend(legend_plot)

# Combine all the trees with the legend
combined_trees_with_legend <- grid.arrange(trees, legend, ncol = 4)
combined_trees_and_maps_with_legend <- grid.arrange(trees, legend, ncol = 4)



# Complete Phylogeny ---------
current_MAPLE_tree <- ape::read.tree("")
current_metadata <- read.csv(file = "")

# Visualizing tree by region
full_tree_tip_region_trait <- data.frame(
  label = current_MAPLE_tree$tip.label,
  region = factor(sequenced_individual_detailed_metadata$region[match(current_MAPLE_tree$tip.label, sequenced_individual_detailed_metadata$strain)]),
  variant = factor(current_metadata$variant[match(current_MAPLE_tree$tip.label, current_metadata$strain)]),
  scorpio_call = factor(current_metadata$scorpio_call[match(current_MAPLE_tree$tip.label, current_metadata$strain)])
)
# Used Taxonium to visualize the full tree
