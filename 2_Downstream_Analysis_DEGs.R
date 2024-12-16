##########################################
# RNA-seq Acylcarnitine Project - Down Stream Analysis (DEGs)
# Author: Kira Liu
# Date: 2024-11-17
# Purpose: Perform downstream analysis on DESeq2 differential expression results, 
#          including filtering, visualization, and exploration of gene overlaps.
# Workflow:
#   1. Load filtered DESeq2 results
#   2. Data inspection and preprocessing
#   3. Filter significant DEGs (differentially expressed genes)
#   4. Analyze overlaps and gene types
#   5. Generate visualizations (e.g., heatmap for selected genes)
##########################################

#### Load Required Libraries ####
library(openxlsx)    # For reading and writing Excel files
library(tidyverse)   # For data manipulation and visualization
library(clusterProfiler)  # For enrichment analysis
library(org.Mm.eg.db)     # For mouse genome annotations
library(pheatmap)    # For heatmap visualization

##########################################
# Step 1: Load DESeq2 Filtered Results
##########################################
data_wd_cold.vs.sc_cold <- read.xlsx("./data/RNAseq_BAT_Filtered_Annotated.xlsx", sheet = 1)
data_wd_cold.vs.wd_rt <- read.xlsx("./data/RNAseq_BAT_Filtered_Annotated.xlsx", sheet = 2)
data_sc_cold.vs.sc_rt <- read.xlsx("./data/RNAseq_BAT_Filtered_Annotated.xlsx", sheet = 3)
data_wd_rt.vs.sc_rt <- read.xlsx("./data/RNAseq_BAT_Filtered_Annotated.xlsx", sheet = 4)
data_interaction_temp_to_diet <- read.xlsx("./data/RNAseq_BAT_Filtered_Annotated.xlsx", sheet = 6)

##########################################
# Step 2: Data Inspection and Preprocessing
##########################################
# Check columns containing '0' values in the first dataset
for (i in colnames(data_wd_cold.vs.sc_cold)) {
  print(paste("Column:", i))
  print(unique(is.na(data_wd_cold.vs.sc_cold[, i])))
  print(unique(data_wd_cold.vs.sc_cold[, i] == 0))
}

# Replace '0' values with NA in the first dataset
data_wd_cold.vs.sc_cold <- data_wd_cold.vs.sc_cold %>%
  mutate(across(everything(), ~ replace(.x, .x == 0, NA)))

##########################################
# Step 3: Define Filtering Function for DEGs
##########################################
# Define thresholds
log2fc_threshold <- 0.6  # Threshold for log2FoldChange
pvalue_threshold <- 0.05  # Threshold for adjusted p-value

# Function to calculate significant DEGs
calculate_degs <- function(data) {
  data %>%
    filter(abs(log2FoldChange) > log2fc_threshold, padj < pvalue_threshold) %>%
    mutate(Significance = case_when(
      log2FoldChange > log2fc_threshold ~ "Upregulated",
      log2FoldChange < -log2fc_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
}

##########################################
# Step 4: Filter DEGs for Each Dataset
##########################################
df_degs_cold <- calculate_degs(data_wd_cold.vs.sc_cold)       # WD_cold vs SC_cold
df_degs_rt <- calculate_degs(data_wd_rt.vs.sc_rt)             # WD_RT vs SC_RT
df_degs_wd <- calculate_degs(data_wd_cold.vs.wd_rt)           # WD_cold vs WD_RT
df_degs_sc <- calculate_degs(data_sc_cold.vs.sc_rt)           # SC_cold vs SC_RT
df_degs_temp_to_diet <- calculate_degs(data_interaction_temp_to_diet)  # Interaction effect

##########################################
# Step 5: Save DEGs to Excel
##########################################
# Create a new Excel workbook and add each DEG dataset
wb <- createWorkbook()
addWorksheet(wb, "DEGs_Cold")
writeData(wb, sheet = "DEGs_Cold", df_degs_cold)
addWorksheet(wb, "DEGs_WD")
writeData(wb, sheet = "DEGs_WD", df_degs_wd)
addWorksheet(wb, "DEGs_SC")
writeData(wb, sheet = "DEGs_SC", df_degs_sc)
addWorksheet(wb, "DEGs_RT")
writeData(wb, sheet = "DEGs_RT", df_degs_rt)
addWorksheet(wb, "DEGs_Temperature_to_Diet")
writeData(wb, sheet = "DEGs_Temperature_to_Diet", df_degs_temp_to_diet)
saveWorkbook(wb, file = "./data/DEGs_Combined.xlsx", overwrite = TRUE)
cat("Excel file has been saved successfully.")

##########################################
# Step 6: Analyze Overlaps and Gene Types
##########################################
# Merge overlapping DEGs between datasets
df_degs_overlaps_diet <- merge(df_degs_cold, df_degs_rt, by = "ENSEMBL")
df_degs_overlaps_temp_to_diet <- merge(df_degs_temp_to_diet, df_degs_overlaps_diet[, c(1, 7, 24)], by = "ENSEMBL")
df_degs_overlaps_temp_to_diet <- df_degs_overlaps_temp_to_diet %>%
  dplyr::select(-matches("\\.x|\\.y"))

# Analyze proportions of GENETYPE
genetype_percentages <- prop.table(table(df_degs_temp_to_diet$GENETYPE)) * 100
genetype_df <- as.data.frame(genetype_percentages)
colnames(genetype_df) <- c("GENETYPE", "Percentage")
print(genetype_df)

##########################################
# Step 7: Heatmap Visualization for Selected Genes
##########################################
# Prepare the expression matrix for heatmap
expression_matrix <- df_degs_overlaps_temp_to_diet[, c(9:28)]
rownames(expression_matrix) <- df_degs_overlaps_temp_to_diet$SYMBOL
expression_matrix <- t(expression_matrix)

# Reorder samples in a specific order
sample_order <- c(
  "WD_cold_1", "WD_cold_2", "WD_cold_3", "WD_cold_4", "WD_cold_5",
  "SC_cold_1", "SC_cold_2", "SC_cold_3", "SC_cold_4", "SC_cold_5",
  "WD_RT_1", "WD_RT_2", "WD_RT_3", "WD_RT_4", "WD_RT_5",
  "SC_RT_1", "SC_RT_2", "SC_RT_3", "SC_RT_4", "SC_RT_5"
)
expression_matrix <- expression_matrix[sample_order, ]

# Plot and save the heatmap
pheatmap(expression_matrix,
         scale = "row",                # Normalize rows
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_rows = FALSE,         # Disable row clustering
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Expression Heatmap for Selected Genes",
         angle_col = 45)               # Rotate column labels

ggsave("./plots/Filtered_DESeq2_Plots(PCA,Heatmaps,MA)/Interaction_Effect/sample_distance_heatmap_Interaction_Effect.png",
       width = 10, height = 8)