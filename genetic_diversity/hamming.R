library(ape)
library(phangorn)
library(readr)

## Use this R script to calculate hamming distances from a folder of fasta files
## Also use later to calculate Tajima's D statistic over time, again from a
## folder of fasta files

setwd(wd)
sequence_directory <- ""
files <- list.files(sequence_directory)
output_directory <- ""
n_files <- length(files)
hamming_output_directory


input_file <- files[1]
input_path <- paste(sequence_directory, input_file, sep = "")
dist_hamming <- as.data.frame(
  as.matrix(
    dist.hamming(
      read.dna(input_path, format = "fasta")
    )
  ),
  ratio = FALSE
)
week_char <- unlist(strsplit(files[1], split = "[._]"))[3]
output_file <- paste(hamming_output_directory, week_char, ".csv", sep = "")
output_path <- paste(output_directory, output_file, sep = "")
write_csv(dist_hamming, output_path)

save_distance_matrix <- function(fasta_file_path, output_path) {
  # Read sequences and calculate Hamming distance matrix using phangorn
  # Save to specified output file

  sequences <- read.dna(fasta_file_path, format = "fasta")
  dist_hamming <- as.data.frame(
    as.matrix(
      dist.hamming(sequences, ratio = FALSE)
    )
  )

  write_csv(dist_hamming, output_path)
}


for (i in 1:n_files) {
  start_time <- Sys.time()
  input_file <- files[i]

  week_char <- unlist(strsplit(input_file, split = "[._]"))[3]


  input_path <- paste(sequence_directory, input_file, sep = "")

  output_file <- paste(hamming_output_directory, week_char, ".csv", sep = "")
  output_path <- paste(output_directory, output_file, sep = "")

  save_distance_matrix(input_path, output_path)
  run_time <- Sys.time() - start_time
  print(paste(
    "Distance matrix", as.character(i), "of", as.character(n_files),
    "done!"
  ))
  print(run_time)
}


run_time <- as.integer((Sys.time() - start_time)[[1]])
print(run_time)






# Hamming Distance matrices
for (i in 1:n_files) {
  start_time <- Sys.time()
  input_file <- files[i]

  week_char <- unlist(strsplit(input_file, split = "[._]"))[3]


  input_path <- paste(sequence_directory, input_file, sep = "")

  sequences <- read.dna(fasta_file_path, format = "fasta")


  output_file <- paste(hamming_output_directory, week_char, ".csv", sep = "")
  output_path <- paste(output_directory, output_file, sep = "")

  save_distance_matrix(input_path, output_path)
  run_time <- Sys.time() - start_time
  print(paste(
    "Distance matrix", as.character(i), "of", as.character(n_files),
    "done!"
  ))
  print(run_time)
}



# Tajima D statistic

library(pegas)
tajima_by_day <- rep(0, n_files)
tajima_pvals_normal <- rep(0, n_files)
tajima_pvals_beta <- rep(0, n_files)

for (i in 1:n_files) {
  start_time <- Sys.time()
  input_file <- files[i]
  day_idx <- strtoi(
    unlist(
      strsplit(input_file, split = "[._]")
    )[3]
  ) + 1 # Plus 1 to move from python indexing to R
  # Read sequences
  input_path <- paste(sequence_directory, input_file, sep = "")
  sequences <- read.dna(input_path, format = "fasta")

  # Calculate and store tajima statistics
  tajima <- tajima.test(sequences)
  tajima_by_day[day_idx] <- tajima$D
  tajima_pvals_normal[day_idx] <- tajima$Pval.normal
  tajima_pvals_beta[day_idx] <- tajima$Pval.beta

  # Print progress
  run_time <- Sys.time() - start_time
  print(paste(
    "Tajima statistic", as.character(i), "of", as.character(n_files),
    "done!"
  ))
  print(run_time)
}


# Save output

tajima_output_path <- ""
day_path <- ""
pvals_normal_path <- ""
pvals_beta_path
write.csv(as.data.frame(tajima_by_day), paste(tajima_output_path, day_path))
write.csv(as.data.frame(tajima_pvals_normal), paste(tajima_output_path, pvals_normal_path))
write.csv(as.data.frame(tajima_pvals_beta), paste(tajima_output_path, pvals_beta_path))
