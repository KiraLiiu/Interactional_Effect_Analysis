##########################################
# Differential Expression Genes (DEGs) and Interaction Effect (IE) Analysis Using DESeq2
# Author: Kira Liu
# Date: 2024-11-18
# Data: Bulk RNA-seq in 4 group
#       Cold Exposure vs Room Temperature (Cold.vs.RT)
#       Western Diet vs Stantard Chow (WD.vs.SC)
# Purpose: To perform differential expression analysis on RNA-seq data
#          focusing on interaction effects between diet and temperature.
# Methods: DESeq2 for differential expression, VST for normalization,
#          ggplot2 for visualization, and AnnotationDbi for gene annotation.
##########################################

# Load required libraries
library(openxlsx)  # For reading Excel files
library(tidyverse) # Data manipulation and visualization
library(DESeq2)    # Differential expression analysis
library(AnnotationDbi)  # Gene annotation
library(org.Mm.eg.db)   # Mouse genome annotation database
library(ggplot2)        # Plotting
library(pheatmap)       # Heatmaps

##########################################
# Step 1: Load raw data
##########################################
# Load RNA-seq count data and sample information
count_data <- read.xlsx("data/raw_data_RNAseq.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE)
sample_info <- read.xlsx("data/raw_data_RNAseq.xlsx", sheet = 3)

# Ensure count data are integers
count_data[1:23] <- sapply(count_data[1:23], as.integer)

##########################################
# Step 2: Prepare metadata and factorize variables
##########################################
# Create combined condition labels
sample_info <- sample_info %>%
  mutate(condition = paste(diet, temperature, sep = "_")) %>%
  group_by(condition) %>%
  mutate(condition = paste(condition, row_number(), sep = "_")) %>%
  ungroup()

# Factorize diet and temperature variables
sample_info$diet <- factor(sample_info$diet, levels = c("SC", "WD"))
sample_info$temperature <- factor(sample_info$temperature, levels = c("RT", "cold"))

# Remove unnecessary samples and rename columns in count data
count_data <- count_data[, -c(4, 16, 20)] # Remove outliers
colnames(count_data) <- sample_info$condition[match(colnames(count_data), sample_info$sample)]

##########################################
# Step 3: Create DESeq2 dataset and perform analysis
##########################################
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ diet + temperature + diet:temperature)

# Filter low count genes
dds <- dds[rowSums(counts(dds)) > 10,]

# Run DESeq2 pipeline
dds <- DESeq(dds)

##########################################
# Step 4: Extract results for specific contrasts
##########################################
# 1. Compare Cold vs RT under WD diet
res_WD <- results(dds, contrast = list(c("temperature_cold_vs_RT", "dietWD.temperaturecold")))

# 2. Compare WD vs SC under Cold exposure
res_cold <- results(dds, contrast = list(c("diet_WD_vs_SC", "dietWD.temperaturecold")))

# 3. Compare WD vs SC under RT
res_RT <- results(dds, contrast = c("diet", "WD", "SC"))

# 4. Compare Cold vs RT under SC diet
res_SC <- results(dds, contrast = c("temperature", "cold", "RT"))

# 5. Assess interaction effect (Temperature x Diet)
res_interaction <- results(dds, name = "dietWD.temperaturecold")

##########################################
# Step 5: Normalize data using VST
##########################################
vsd_transformed <- vst(dds, blind = FALSE)
normalized_counts <- as.data.frame(assay(vsd_transformed))
normalized_counts$ENSEMBL <- rownames(normalized_counts) # Add ENSEMBL IDs

##########################################
# Step 6: Gene annotation
##########################################
# Annotate genes using ENSEMBL IDs
res_interaction_df <- as.data.frame(res_interaction)
res_interaction_df$ENSEMBL <- rownames(res_interaction_df)
gene_annotations <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = res_interaction_df$ENSEMBL,
                                          columns = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                                          keytype = "ENSEMBL",
                                          multiVals = "first")

# Merge annotations and normalized data
annotated_res_interaction_df <- res_interaction_df %>%
  merge(gene_annotations, by = "ENSEMBL", all = TRUE) %>%
  merge(normalized_counts, by = "ENSEMBL", all = TRUE)

# Export annotated results
write.xlsx(annotated_res_interaction_df, file = "data/DESeq2_Interaction_Effect_results_annotated.xlsx", overwrite = TRUE)

##########################################
# Step 7: Visualization
##########################################
# MA plots for contrasts
png("./plots/MA_Plot_WD_vs_SC_in_cold.png", width = 800, height = 600)
plotMA(res_cold, main = "WD vs SC in Cold Environment", ylim = c(-5, 5))
dev.off()

png("./plots/MA_Plot_WD_vs_SC_in_RT.png", width = 800, height = 600)
plotMA(res_RT, main = "WD vs SC in RT Environment", ylim = c(-5, 5))
dev.off()

png("./plots/MA_Plot_Interaction_Effect.png", width = 800, height = 600)
plotMA(res_interaction, main = "Interaction Effect (Temperature x Diet)", ylim = c(-5, 5))
dev.off()

# Sample distance heatmap
sample_distances <- dist(t(assay(vsd_transformed)))
png("./plots/sample_distance_heatmap_Interaction_Effect.png", width = 800, height = 800)
pheatmap(as.matrix(sample_distances), clustering_distance_rows = sample_distances, 
         clustering_distance_cols = sample_distances, main = "Sample Distance Heatmap")
dev.off()

# PCA plot for sample clustering
png("./plots/PCA_plot_Interaction_Effect.png", width = 800, height = 800)
plotPCA(vsd_transformed, intgroup = c("diet", "temperature")) + ggtitle("PCA Plot")
dev.off()