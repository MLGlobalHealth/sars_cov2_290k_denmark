library("dplyr")
library("lubridate")
library("synthpop")

# Read the true data
df <- read.csv("")

# Add a days column for synthetic date generation
df$days <- as.numeric(
  difftime(as.Date(df$date), as.Date("2021-01-01"), units = "days")
)

sub_df <- subset(df, select = -c(Unnamed..0, PERSON_ID, date, age_groups, branch_lengths, month))

# Density smoothing + increased ctree.minbucket decrease disclosure risk
syn_data <- syn(sub_df, smoothing = "density", ctree.minbucket = 10)

# Safety check
replicated.uniques(syn_data, sub_df)
sdc(syn_data, sub_df)

# Quality check
compare(syn_data, sub_df, stat = "counts")

# Augment
syn_df <- syn_data$syn

## Add branch lengths
syn_df$syn$branch_lengths <- syn_df$syn$rate * syn_df$syn$branch_lengths_time

## Add PERSON ID
set.seed(123) # For reproducibility
syn_df$PERSON_ID <- sample(100000:999999, size = nrow(syn_df), replace = FALSE)

## Add date
syn_df$date <- start_date + days(syn_df$days)

syn_df$month <- month(syn_df$date)

### Drop the 'days' column
syn_df <- syn_df %>% select(-days)

## Add age groups
### Define age groups
age_groups <- c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")

syn_df <- syn_df %>% mutate(age_groups = age_groups[pmin((ages %/% 15) + 1, 6)])

# Save
write.csv(syn_df, "data/synthetic/maple_tip_lengths.csv")
