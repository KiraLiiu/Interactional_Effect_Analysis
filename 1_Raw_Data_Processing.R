##########################################
# RNA-seq Acylcarnitine Project - Row Data Processing (DESeq2)
# Author: Kira Liu
# Date: 2024-11-17
# Purpose: Perform RNA-seq differential expression analysis for specific 
#          comparisons between diet and temperature conditions.
# Workflow:
#   1. Load RNA-seq count data and metadata
#   2. Preprocess and filter data
#   3. Define and run modularized DESeq2 analysis for condition comparisons
#   4. Save results for further downstream analysis
##########################################

# Load required libraries
library(openxlsx)  # For reading Excel files
library(tidyverse) # Data manipulation and visualization
library(dplyr)     # Enhanced data manipulation

##########################################
# Step 1: Load raw data
##########################################
# Load RNA-seq count data (rows: genes, columns: samples)
count_data <- read.xlsx("data/raw_data_RNAseq.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE)

# Load sample information (metadata for samples)
sample_info <- read.xlsx("data/raw_data_RNAseq.xlsx", sheet = 3)

##########################################
# Step 2: Preprocess count data
##########################################
# Ensure that all count data columns are integers (required by DESeq2)
count_data[1:23] <- sapply(count_data[1:23], as.integer)

##########################################
# Step 3: Add condition labels
##########################################
# Create a new column in the sample metadata combining diet and temperature
sample_info <- sample_info %>%
  mutate(condition = paste(diet, temperature, sep = "_"))

##########################################
# Step 4: Function to perform DESeq2 analysis for specific comparisons
##########################################
run_deseq2_comparison <- function(count_data, sample_info, condition_filter, condition_levels, result_name) {
  # Purpose: Perform DESeq2 analysis for specific conditions
  # Arguments:
  #   - count_data: Raw RNA-seq count matrix
  #   - sample_info: Metadata with sample conditions
  #   - condition_filter: Filter for specific conditions (e.g., "cold")
  #   - condition_levels: Factor levels for comparison (e.g., c("SC_cold", "WD_cold"))
  #   - result_name: Variable name to store the results in the global environment
  
  # Step 1: Filter metadata based on the specified condition
  selected_samples <- sample_info %>%
    filter(str_detect(condition, condition_filter)) %>% # Filter rows based on condition
    dplyr::select(sample, condition)                   # Select relevant columns
  
  # Step 2: Subset count data for the selected samples
  filtered_counts <- count_data %>%
    dplyr::select(all_of(selected_samples$sample))      # Select columns matching sample IDs
  
  # Step 3: Run DESeq2 analysis
  # Call the DESeq2 wrapper function for DEA
  result <- run_DESeq2_analysis(
    count_data = filtered_counts,
    sample_info = selected_samples,
    condition_levels = condition_levels
  )
  
  # Step 4: Save results in the global environment for later use
  assign(result_name, result, envir = .GlobalEnv)
}

##########################################
# Step 5: Perform DESeq2 analyses for specific condition comparisons
##########################################

# 1. Compare WD_cold vs SC_cold
run_deseq2_comparison(
  count_data = count_data,
  sample_info = sample_info,
  condition_filter = "cold",
  condition_levels = c("SC_cold", "WD_cold"),
  result_name = "data_bat_rna_cold"
)

# 2. Compare WD_RT vs SC_RT
run_deseq2_comparison(
  count_data = count_data,
  sample_info = sample_info,
  condition_filter = "RT",
  condition_levels = c("SC_RT", "WD_RT"),
  result_name = "data_bat_rna_RT"
)

# 3. Compare WD_cold vs WD_RT
run_deseq2_comparison(
  count_data = count_data,
  sample_info = sample_info,
  condition_filter = "WD",
  condition_levels = c("WD_RT", "WD_cold"),
  result_name = "data_bat_rna_WD"
)

# 4. Compare SC_cold vs SC_RT
run_deseq2_comparison(
  count_data = count_data,
  sample_info = sample_info,
  condition_filter = "SC",
  condition_levels = c("SC_RT", "SC_cold"),
  result_name = "data_bat_rna_SC"
)

##########################################
# Step 6: Special case comparison (WD_cold vs SC_RT)
##########################################

# Manually filter for WD_cold and SC_RT samples
selected_samples <- sample_info[6:15, ][, c(1, 4)] # Subset specific rows and columns

# Subset count data for selected samples
filtered_counts <- count_data %>%
  dplyr::select(all_of(selected_samples$sample))      # Match columns to sample IDs

# Run DESeq2 analysis directly for this special case
data_bat_rna <- run_DESeq2_analysis(
  count_data = filtered_counts,
  sample_info = selected_samples,
  condition_levels = c("SC_RT", "WD_cold")
)

##########################################
# Step 7: Clean up intermediate variables
##########################################
# Remove all variables from the current environment (optional)
rm(list = ls())