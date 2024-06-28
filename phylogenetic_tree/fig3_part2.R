# Analysis of demographic details and ascertainment bias

# Author: Mark Khurana (mark.khurana@sund.ku.dk)

# Importing packages
library(
  ape, phytools, TreeTools, dplyr, tidyverse,
  data.table, dbplyr, lubridate, rlang, foreach,
  doParallel, DSTora, ROracle, DSTcolectica, DSTdb,
  DBI, parallel, ggsignif, Rcpp
)

# Establish connection
drv <- dbDriver("Oracle")
conn <- DSTora::OraGenvej("", dbuser = "")
con2 <- DSTora::OraGenvej("", dbuser = "")

# To access data tables
dbListTables(conn)

covid_test <- dbReadTable()

lifelines <- dbReadTable()

lifelines_koen <- dbReadTable()

query <- ""
geo_address <- dbGetQuery(conn = con2, statement = query)

general_address <- dbReadTable()


# Data Preparation


# All sequenced individuals, regardless of whether
# we have full metadata for them or not
all_sequenced_individuals <- read.csv(file = "")
temp_data <- left_join(all_sequenced_individuals,
  lifelines %>% distinct(PERSON_ID, .keep_all = TRUE),
  by = "PERSON_ID"
)

all_sequenced_individuals_combined <- left_join(temp_data,
  lifelines_koen %>% distinct(PERSON_ID, .keep_all = TRUE),
  by = "PERSON_ID"
)

all_sequenced_individuals_sex <- left_join(all_sequenced_individuals,
  lifelines_koen %>% distinct(PERSON_ID, .keep_all = TRUE),
  by = "PERSON_ID"
)

all_sequenced_individuals_combined <- all_sequenced_individuals_combined %>%
  mutate(
    BIRTHDAY = ymd(BIRTHDAY),
    date = ymd(date),
    age_at_infection = as.numeric(
      difftime(date, BIRTHDAY, units = "days")
    ) / 365.25
  )

all_sequenced_individuals_combined_unique <- all_sequenced_individuals_combined %>%
  mutate(
    BIRTHDAY = ymd(BIRTHDAY),
    date = ymd(date),
    age_at_infection = as.numeric(difftime(date, BIRTHDAY, units = "days")) / 365.25
  ) %>%
  filter(!is.na(PERSON_ID)) %>%
  distinct(PERSON_ID, .keep_all = TRUE)

# 293287 individuals in total
length(all_sequenced_individuals$PERSON_ID[!is.na(all_sequenced_individuals$PERSON_ID)])
# 292481 unique individuals
length(all_sequenced_individuals_combined_unique$PERSON_ID[!is.na(all_sequenced_individuals$PERSON_ID)])

# Adding geographic data
sequenced_individual_detailed_metadata <- all_sequenced_individuals_combined_unique %>%
  left_join(geo_address, by = "PERSON_ID")
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  left_join(general_address, by = "ID")
# Filter rows based on date range
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  mutate(
    BOP_VFRA = as.Date(BOP_VFRA),
    BOP_VTIL = as.Date(BOP_VTIL)
  )
# Filter to ensure that individuals were living at the residence during their positive test
sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata %>%
  filter(date >= BOP_VFRA & date <= BOP_VTIL) %>%
  distinct(PERSON_ID, .keep_all = TRUE)

# Metadata for all individuals with a positive COVID test in the given period

all_positive_individuals <- covid_test %>%
  filter(
    !is.na(PERSON_ID),
    between(PRDATE, as.Date("2021-01-01"), as.Date("2021-12-31"))
  ) %>%
  filter(SVARRESULTAT == 1)

# Counting number of positive tests
length(all_positive_individuals$PERSON_ID) # 966094 positive tests

all_positive_pcr_individuals_missing_id <- all_positive_individuals %>%
  filter(CASEDEF == "SARS2") %>%
  filter(SVARRESULTAT == 1)
length(all_positive_pcr_individuals_missing_id$PERSON_ID) # 731122 positive PCR tests

# Creating metadata file
common_person_ids <- all_positive_individuals %>%
  select(PERSON_ID)
# Step 1: Filter 'lifelines' for rows where 'PERSON_ID' is in 'all_positive_individuals'
filtered_lifelines <- lifelines %>%
  semi_join(common_person_ids, by = "PERSON_ID")
# Step 2: Join 'filtered_lifelines' with 'lifelines' on 'PERSON_ID' and take the first 'BIRTHDAY'
merged_data <- all_positive_individuals %>%
  left_join(
    filtered_lifelines %>%
      group_by(PERSON_ID) %>%
      summarize(BIRTHDAY = first(BIRTHDAY)),
    by = "PERSON_ID"
  )
# Step 3: Join 'merged_data' with 'lifelines_koen' on 'PERSON_ID'
new_merged_metadata <- merged_data %>%
  left_join(
    lifelines_koen %>%
      group_by(PERSON_ID) %>%
      summarize(KOEN = first(KOEN)),
    by = "PERSON_ID"
  )

positive_individual_basic_metadata <- new_merged_metadata
positive_individual_basic_metadata <- positive_individual_basic_metadata %>%
  mutate(age_at_infection = as.numeric(difftime(PRDATE, BIRTHDAY, units = "days")) / 365.25) %>%
  distinct(PERSON_ID, .keep_all = TRUE)

# Step 4: Join data with 'geo_address' on 'PERSON_ID'
positive_individual_detailed_metadata <- new_merged_metadata %>%
  left_join(geo_address, by = "PERSON_ID")
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  left_join(general_address, by = "ID")
# Step 5: Filter rows based on date range
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  mutate(
    BOP_VFRA = as.Date(BOP_VFRA),
    BOP_VTIL = as.Date(BOP_VTIL),
    BIRTHDAY = as.Date(BIRTHDAY),
    PRDATE = as.Date(PRDATE)
  )
# Proceed with the filter
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  filter(PRDATE >= BOP_VFRA & PRDATE <= BOP_VTIL)

# Step 6: Only keep one row per person:
# Arrange the data by 'PRDATE' in ascending order within each 'PERSON_ID' group
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  arrange(PERSON_ID, PRDATE)
# Keep only the rows with the first/lowest 'PRDATE'
# within each 'PERSON_ID' group - i.e, only unique individuals
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  distinct(PERSON_ID, .keep_all = TRUE)

# Step 7: Add an age_at_infection column
positive_individual_detailed_metadata <- positive_individual_detailed_metadata %>%
  mutate(age_at_infection = as.numeric(difftime(PRDATE, BIRTHDAY, units = "days")) / 365.25)

# Characteristics of sequenced vs non-sequenced individuals -------

# Age ----
length(positive_individual_basic_metadata$age_at_infection) # 686,851 unique individuals
sum(is.na(positive_individual_basic_metadata$age_at_infection)) # 1798 missing
age_groups <- cut(positive_individual_basic_metadata$age_at_infection,
  breaks = c(0, 15, 30, 45, 60, 75, Inf),
  labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"),
  include.lowest = TRUE
)
positive_df <- data.frame(age_groups)
table_positive <- table(positive_df$age_groups)
percentage_positive <- prop.table(table_positive) * 100
print(table_positive)
print(percentage_positive)
# Repeat the process for sequenced_individuals
length(all_sequenced_individuals_combined_unique$age_at_infection) # 292,481 unique individuals
sum(is.na(all_sequenced_individuals_combined_unique$age_at_infection)) # 898 missing
age_groups_sequenced <- cut(all_sequenced_individuals_combined_unique$age_at_infection,
  breaks = c(0, 15, 30, 45, 60, 75, Inf),
  labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"),
  include.lowest = TRUE
)
sequenced_df <- data.frame(age_groups_sequenced)
table_sequenced <- table(sequenced_df$age_groups_sequenced)
percentage_sequenced <- prop.table(table_sequenced) * 100
# Display table for sequenced_individuals
print(table_sequenced)
print(percentage_sequenced)
# Perform statistical test
chi_square_result <- chisq.test(table_positive, table_sequenced)
print(chi_square_result)


# Sex ----
positive_individual_basic_metadata$KOEN <- as.factor(positive_individual_basic_metadata$KOEN)
all_sequenced_individuals_combined_unique$KOEN <- as.factor(
  all_sequenced_individuals_combined_unique$KOEN
)

summary(positive_individual_basic_metadata$KOEN)
summary(all_sequenced_individuals_combined_unique$KOEN)
positive_counts <- table(positive_individual_basic_metadata$KOEN)
sequenced_counts <- table(all_sequenced_individuals_combined_unique$KOEN)
prop.table(positive_counts) * 100
prop.table(sequenced_counts) * 100
prop_test_result <- prop.test(positive_counts, sequenced_counts)
print(prop_test_result)


# Region ----
positive_individual_detailed_metadata$REGIONSKODE <- as.factor(
  positive_individual_detailed_metadata$REGIONSKODE
)
sequenced_individual_detailed_metadata$REGIONSKODE <- as.factor(
  sequenced_individual_detailed_metadata$REGIONSKODE
)
# Breaking down by region:
positive_counts <- table(positive_individual_detailed_metadata$REGIONSKODE)
sequenced_counts <- table(sequenced_individual_detailed_metadata$REGIONSKODE)
positive_percentage <- prop.table(positive_counts) * 100
sequenced_percentage <- prop.table(sequenced_counts) * 100
contingency_table <- matrix(c(positive_counts, sequenced_counts), ncol = 2)
# Perform chi-square test
chi_square_result <- chisq.test(contingency_table)


# Figure with Demographic Breakdown, Part of Figure 3 -------

# Age Distribution
all_sequenced_individuals_combined$age_group <- cut(
  all_sequenced_individuals_combined$age_at_infection,
  breaks = c(0, 15, 30, 45, 60, 75, Inf),
  labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"),
  include.lowest = TRUE
)

all_sequenced_individuals_combined$age_group <- factor(
  all_sequenced_individuals_combined$age_group,
  levels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")
)
all_sequenced_individuals_combined <- all_sequenced_individuals_combined[
  !is.na(all_sequenced_individuals_combined$age_group),
]

age_plot <- ggplot(all_sequenced_individuals_combined, aes(x = age_group, fill = age_group)) +
  geom_bar(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "plasma") + # Use plasma viridis colors for fill
  labs(x = "Age Group", y = "Count") +
  theme_classic()

# Region
region_name_mapping <- c(
  "1081" = "N",
  "1082" = "M",
  "1083" = "SY",
  "1084" = "H",
  "1085" = "SJ"
)
sequenced_individual_detailed_metadata$region <- factor(
  region_name_mapping[as.character(sequenced_individual_detailed_metadata$REGIONSKODE)]
)

sequenced_individual_detailed_metadata <- sequenced_individual_detailed_metadata[
  !is.na(sequenced_individual_detailed_metadata$region),
]

region_plot <- ggplot(
  sequenced_individual_detailed_metadata,
  aes(x = region, fill = region)
) +
  geom_bar(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) + # Use viridis colors for fill
  labs(x = "Region", y = "Count") +
  theme_classic()

# Major Variant
sequenced_individuals <- sequenced_individuals[!is.na(sequenced_individuals$variant), ]

sequenced_individuals <- sequenced_individuals %>%
  mutate(major_variant = ifelse(grepl("^Omicron", variant), "Omicron", variant))

major_variant_plot <- ggplot(
  sequenced_individuals,
  aes(x = major_variant, fill = major_variant)
) +
  geom_bar(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  labs(x = "Major Variant", y = "Count") +
  theme_classic()

demographics_plot <- grid.arrange(major_variant_plot, age_plot, region_plot,
  ncol = 1
)

# Save the combined figure with legend
# ggsave("", demographics_plot, width = 4, height = 9)
