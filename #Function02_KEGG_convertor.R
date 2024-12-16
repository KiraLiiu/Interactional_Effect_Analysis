##########################################
# Function: Convert ENTREZIDs to Gene Names
# Author: Kira Liu
# Date: 2024-10-29
# Purpose: To map ENTREZIDs in KEGG pathway analysis results to corresponding gene names.
# Inputs:
#   - kegg_result: Data frame containing KEGG pathway analysis results.
#   - mapping_df: Data frame with ENTREZID-to-Gene.Name mappings.
#   - entrez_col: Column name in kegg_result containing ENTREZIDs (default = "geneID").
#   - id_col: Column name in mapping_df containing ENTREZIDs (default = "ENTREZID").
#   - name_col: Column name in mapping_df containing gene names (default = "Gene.Name").
# Outputs:
#   - Updated KEGG result data frame with ENTREZIDs converted to gene names.
##########################################

convert_entrez_to_name <- function(kegg_result, mapping_df, entrez_col = "geneID", id_col = "ENTREZID", name_col = "Gene.Name") {
  
  # Step 1: Check if required columns exist in the mapping data frame
  if (!all(c(id_col, name_col) %in% colnames(mapping_df))) {
    stop("The mapping data frame must contain the specified ID and Name columns.")
  }
  
  # Step 2: Apply mapping to the specified ENTREZID column in the KEGG result
  kegg_result[[entrez_col]] <- sapply(kegg_result[[entrez_col]], function(entrez_ids) {
    # Split the ENTREZID string into individual IDs
    ids <- unlist(strsplit(entrez_ids, "/"))
    
    # Map ENTREZIDs to Gene Names using the provided mapping data frame
    gene_names <- mapping_df %>%
      filter(!!sym(id_col) %in% ids) %>%  # Dynamically select the ID column
      pull(!!sym(name_col))               # Dynamically select the Name column
    
    # Combine the gene names into a single string, separated by "/"
    paste(gene_names, collapse = "/")
  })
  
  # Step 3: Optional - Clean up the "Description" column if it exists
  if ("Description" %in% colnames(kegg_result)) {
    # Remove redundant information from the Description column
    kegg_result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", 
                                    kegg_result$Description, perl = TRUE)
  }
  
  # Step 4: Return the updated KEGG result data frame
  return(kegg_result)
}

# Example Usage:
# Assuming `kegg_results_df` contains KEGG results and `gene_mapping_df` provides ENTREZID-to-Gene.Name mapping:
# updated_kegg_results <- convert_entrez_to_name(kegg_results_df, gene_mapping_df)