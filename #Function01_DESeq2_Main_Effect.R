##########################################
# RNA-seq Differential Expression Analysis Function Using DESeq2
# Author: Kira Liu
# Date: 2024-11-13
# Purpose: To automate RNA-seq differential expression analysis and visualization 
#          using DESeq2, with support for gene annotation and customized outputs.
# Methods: DESeq2 for differential expression analysis, VST for normalization,
#          AnnotationDbi for gene annotation, and ggplot2/pheatmap for visualization.
# Workflow: 
#   1. Construct DESeq2 dataset
#   2. Filter low-count genes
#   3. Run DESeq2 pipeline
#   4. Annotate results with gene information
#   5. Normalize data using VST
#   6. Generate output files (Excel and plots)
##########################################

run_DESeq2_analysis <- function(count_data, sample_info, condition_levels = c("SC_cold", "WD_cold"), alpha_level = 0.05) {
  
  # Load required libraries
  library(DESeq2)        # For differential expression analysis
  library(AnnotationDbi) # For gene annotation
  library(org.Mm.eg.db)  # Mouse genome annotation database
  library(pheatmap)      # For heatmap visualization
  library(ggplot2)       # For PCA and other plots
  
  ##########################################
  # Step 1: Prepare DESeq2 dataset
  ##########################################
  # Create a DESeq2 dataset using count data and sample metadata
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_info,
                                design = ~ condition)
  
  ##########################################
  # Step 2: Filter low-count genes
  ##########################################
  # Remove genes with low counts (sum of counts > 10 across samples)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  # Ensure condition levels are correctly ordered
  dds$condition <- factor(dds$condition, levels = condition_levels)
  
  ##########################################
  # Step 3: Run DESeq2 differential expression analysis
  ##########################################
  # Perform DESeq2 analysis on the filtered dataset
  dds <- DESeq(dds)
  
  ##########################################
  # Step 4: Extract and format results
  ##########################################
  # Retrieve the names of results (contrasts)
  results_names <- resultsNames(dds)
  
  # Extract results for the specified condition comparison
  res <- results(dds, name = results_names[2], alpha = alpha_level)
  
  # Convert results to a dataframe and add ENSEMBL IDs
  res_df <- as.data.frame(res)
  res_df$ENSEMBL <- rownames(res_df)
  
  ##########################################
  # Step 5: Normalize data using variance stabilizing transformation (VST)
  ##########################################
  # Normalize count data for better visualization and comparison
  vsd <- vst(dds, blind = FALSE)
  normalized_data <- as.data.frame(assay(vsd))
  
  # Add sample condition names to column headers
  colnames(normalized_data) <- sample_info$condition[match(colnames(normalized_data), sample_info$sample)]
  normalized_data$ENSEMBL <- rownames(normalized_data)
  
  ##########################################
  # Step 6: Annotate results with gene information
  ##########################################
  # Retrieve gene annotations using ENSEMBL IDs
  annotations <- AnnotationDbi::select(org.Mm.eg.db,
                                       keys = res_df$ENSEMBL,
                                       columns = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                                       keytype = "ENSEMBL",
                                       multiVals = "first")
  
  # Merge results, annotations, and normalized data
  annotated_results <- merge(res_df, annotations, by = "ENSEMBL", all = TRUE)
  full_data <- merge(annotated_results, normalized_data, by = "ENSEMBL", all = TRUE)
  
  ##########################################
  # Step 7: Export results to an Excel file
  ##########################################
  # Save annotated results to an Excel file
  output_file <- paste0("./data/RNAseq_", results_names[2], ".xlsx")
  write.xlsx(full_data, file = output_file, overwrite = TRUE)
  
  ##########################################
  # Step 8: Generate MA plot
  ##########################################
  # Create an MA plot to visualize log fold changes vs mean expression
  ma_plot_file <- paste0("./plots/MA_plot_", results_names[2], ".png")
  png(ma_plot_file)
  plotMA(res, ylim = c(-5, 5), main = "MA Plot")
  dev.off()
  
  ##########################################
  # Step 9: Generate sample distance heatmap
  ##########################################
  # Calculate Euclidean distances between samples based on normalized counts
  sample_dist <- dist(t(assay(vsd)))
  heatmap_file <- paste0("./plots/sample_distance_heatmap_", results_names[2], ".png")
  png(heatmap_file)
  pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, 
           clustering_distance_cols = sample_dist, main = "Sample Distance Heatmap")
  dev.off()
  
  ##########################################
  # Step 10: Generate PCA plot
  ##########################################
  # Perform PCA and create a scatter plot for sample clustering
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_plot_file <- paste0("./plots/PCA_plot_", results_names[2], ".png")
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Customize and save PCA plot
  ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(pca_plot_file)
  
  ##########################################
  # Return full annotated results for further analysis
  ##########################################
  return(full_data)
}

# Example usage of the function
# Assuming `rawdata_bat_rna_cold` is the count data matrix
# and `sample_info_cold` is the sample metadata
data_bat_rna_cold <- run_DESeq2_analysis(
  count_data = rawdata_bat_rna_cold,
  sample_info = sample_info_cold,
  condition_levels = c("SC_cold", "WD_cold")
)