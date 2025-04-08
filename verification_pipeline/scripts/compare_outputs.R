#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Function to load data (no header, tab-delimited)
load_data <- function(file, skip_header = FALSE) {
  if (skip_header) {
    return(read.table(file, header=FALSE, sep="\t", skip=1))
  } else {
    return(read.table(file, header=FALSE, sep="\t"))
  }
}

# Function to compare outputs
compare_outputs <- function(file1, file2, fig_output) {
  
  # Check if the files are identical
  # if (identical(readLines(file1), readLines(file2))) {
  #   cat("\nThe files are identical.\n")
  #   return()  # Skip further processing if files are identical
  # }

  # Determine if the files are "ihh12" output and remove the header if needed
  skip_header1 <- grepl("ihh12", file1)
  skip_header2 <- grepl("ihh12", file2)
  
  # Load the two files (no header), skip the header for ihh12 files
  df1 <- load_data(file1, skip_header = skip_header1)
  df2 <- load_data(file2, skip_header = skip_header2)
  
  # Ensure the first column is used as the key for merging (use the first column as 'CHR')
  colnames(df1) <- c("CHR", paste0("V", 2:ncol(df1)))
  colnames(df2) <- c("CHR", paste0("V", 2:ncol(df2)))
  
  # Merge the data by the first column ('CHR')
  merged <- merge(df1, df2, by="CHR", suffixes=c("_sc", "_sb"))
  
  # Compare the last columns (selection statistic) in both data frames
  last_col_sc <- ncol(df1)
  last_col_sb <- ncol(df2)
  
  # Calculate the difference between the selection statistics
  merged$diff <- merged[, last_col_sc] - merged[, last_col_sb]
  
  # Print the first few rows of the result (CHR and diff)
  print(head(merged[, c("CHR", "diff")]))
  
  # Print summary statistics for the differences
  print(summary(merged$diff))
  
  # Calculate Pearson correlation between the selection statistics
  correlation <- cor(merged[, last_col_sc], merged[, last_col_sb], method="pearson")
  cat("\nPearson correlation between the two selection statistics: ", correlation, "\n")

  # Calculate Spearman's rank correlation
spearman_correlation <- cor(merged[, last_col_sc], merged[, last_col_sb], method="spearman")
cat("\nSpearman's rank correlation: ", spearman_correlation, "\n")

# Calculate Kendall's Tau correlation
kendall_correlation <- cor(merged[, last_col_sc], merged[, last_col_sb], method="kendall")
cat("\nKendall's Tau correlation: ", kendall_correlation, "\n")
  
  # Percentage of common values retained in `sc`
  common_chr <- intersect(df1$CHR, df2$CHR)
  total_chr_sc <- nrow(df1)
  total_chr_sb <- nrow(df2)
  
  common_percentage_sc <- length(common_chr) / total_chr_sc * 100
  common_percentage_sb <- length(common_chr) / total_chr_sb * 100
  
  # Percentage of uncommon values excluded from both `sc` and `sb`
  uncommon_chr_sc <- setdiff(df1$CHR, df2$CHR)
  uncommon_chr_sb <- setdiff(df2$CHR, df1$CHR)
  
  excluded_percentage_sc <- length(uncommon_chr_sc) / total_chr_sc * 100
  excluded_percentage_sb <- length(uncommon_chr_sb) / total_chr_sb * 100
  
  # Print the percentages
  cat("\nPercentage of common values retained in SC: ", common_percentage_sc, "%\n")
  cat("Percentage of common values retained in SB: ", common_percentage_sb, "%\n")
  cat("Percentage of uncommon values excluded from SC: ", excluded_percentage_sc, "%\n")
  cat("Percentage of uncommon values excluded from SB: ", excluded_percentage_sb, "%\n")
  
  # Plot the correlation
  plot <- ggplot(merged, aes(x = merged[, last_col_sc], y = merged[, last_col_sb])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "Correlation between SC and SB selection statistics",
         x = "SC Statistic", y = "SB Statistic") +
    theme_minimal()
  
  # Save the plot
  ggsave(fig_output, plot)
  cat("\nCorrelation plot saved as", fig_output, "\n")
}

# Main execution
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) {
  stop("Usage: compare_outputs.R <run1_output> <run2_output> <output_figure>")
}

compare_outputs(args[1], args[2], args[3])
