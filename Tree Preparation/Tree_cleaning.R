# Code to identify outliers

library(
  "ape", "adegenet", "TreeTools", "tidyverse", "phangorn", "patchwork",
  "phytools", "phylotools", "ggplot2", "seqinr", "readr", "purrr", "openxlsx",
  "lubridate", "dplyr", "data.table", "magrittr", "tidyr", "parallel",
  "reshape2", "pegas", "stats", "lubridate", "mlesky", "MASS", "goft",
  "fitdistrplus", "brms"
)

# Reading metadata and matching names to sequences in order to get dates
summary_data <- read.csv()
metadata <- read.csv()
colnames(summary_data)[colnames(summary_data) == "file_name"] <- "ConsensusFilename"
metadata$ConsensusFilename <- gsub("\\.fa$|\\.fasta$", "", metadata$ConsensusFilename)
merged_metadata <- merge(
  summary_data,
  metadata,
  by = "ConsensusFilename",
  all.x = TRUE,
  all.y = TRUE
)
merged_metadata$seq_name <- gsub("\\.fa$|\\.fasta$", "", merged_metadata$seq_name)
merged_metadata$strain <- merged_metadata$seq_name

# Importing pango lineages and adding to merged_metadata
consensus_pango <- read.csv("")
colnames(consensus_pango)[1] <- "strain"
merged_metadata <- merge(
  merged_metadata,
  consensus_pango,
  by = "strain",
  all.x = TRUE,
  all.y = TRUE
)


# MAPLE run for all sequences
# Finding outliers for branch lengths and removing the reference sequence ---------
# Readin in MAPLE tree
tr2 <- read.tree("newids_consensus_2021_tree.tree")

# Fit gamma distribution and then prune outliers
# Extract branch lengths from the tree
branch_lengths <- as.vector(tr2$edge.length)
branch_lengths <- branch_lengths[is.finite(branch_lengths)]
# Remove zero branch lengths
non_zero_branch_lengths <- data.frame(branch_lengths[branch_lengths > 0])
nonzerobl <- branch_lengths[branch_lengths > 0]
fit_gamma <- fitdistr(nonzerobl, "gamma")
# Determine a cutoff based on the distribution (e.g., 99th percentile)
cutoff_percentile <- qgamma(0.99, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2])
# Prune branches exceeding the cutoff
tr2$tip.label[tr2$edge.length > cutoff_percentile]
pruned_tree <- drop.tip(tr2, tr2$tip.label[tr2$edge.length > cutoff_percentile])
ape::write.tree(pruned_tree, file = "", digits = 16)


# BRMS --------
# Read in chronumental tree and output metadata
tr_chron <- read.tree("")
tsv_data <- read.table("", header = TRUE, sep = "\t")
chronumental_date_predicted_merged_data <- merge(merged_metadata, tsv_data, by = "strain")
chronumental_date_predicted_merged_data$predicted_date <- decimal_date(
  as.Date(chronumental_date_predicted_merged_data$predicted_date)
)
chronumental_date_predicted_merged_data$date <- decimal_date(
  as.Date(chronumental_date_predicted_merged_data$date)
)

# Running BRMS (Bayesian Regression Models using Stan) to see which to remove as outliers
data <- chronumental_date_predicted_merged_data[, c("strain", "date", "predicted_date", "variant")]

library(rstanarm)
d <- data

m <- stan_glm(formula = predicted_date ~ -1 + date, data = d, adapt_delta = 0.99)
draws <- as.data.frame(m)
coeff <- draws[, 1]
low <- rep(0, nrow(draws))
high <- rep(0, nrow(draws))
for (i in seq_len(nrow(d))) {
  low[i] <- quantile(d$date[i] * coeff, probs = 0.025)
  high[i] <- quantile(d$date[i] * coeff, probs = 0.975)
}
out <- d$predicted_date > low & d$predicted_date > high

100 * sum(out) / nrow(d)

sequences_to_prune <- d$strain[out] # Extract names of sequences outside the interval
pruned_phy_tree <- drop.tip(tr_chron, sequences_to_prune)
ape::write.tree(pruned_phy_tree, file = "", digits = 16)
