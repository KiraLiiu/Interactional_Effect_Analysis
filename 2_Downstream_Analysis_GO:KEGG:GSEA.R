##########################################
# RNA-seq Acylcarnitine Project - Down Stream Analysis (GO/KEGG/GSEA/TRNs)
# Author: Kira Liu
# Date: 2020-11-18
# Purpose: Perform GO and KEGG enrichment analysis, GSEA, and transcription factor (TF) network construction
# Workflow:
#   1. Enrichment analysis (GO and KEGG) for significant DEGs
#   2. Gene Set Enrichment Analysis (GSEA) for GO/KEGG terms
#   3. Extract gene lists for transcription factor (TF) analysis
#   4. Construct and export TF-Gene regulatory networks
##########################################

#### Load Required Libraries ####
library(clusterProfiler) # For GO/KEGG and GSEA analysis
library(org.Mm.eg.db)    # Mouse genome annotations
library(dplyr)           # Data manipulation
library(tidyr)           # For data reshaping
library(openxlsx)        # For reading and writing Excel files
library(tidygraph)       # For network construction
library(ggraph)          # For network visualization

##########################################
# Step 1: GO Enrichment Analysis for DEGs
##########################################

# Define a function for GO enrichment
perform_go_enrichment <- function(degs, universe, orgdb, ont = "BP", q_cutoff = 0.05, p_cutoff = 0.05) {
  degs %>%
    clusterProfiler::enrichGO(
      OrgDb = orgdb,
      ont = ont,
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      universe = universe,
      readable = TRUE,
      pAdjustMethod = "BH"
    ) %>%
    clusterProfiler::simplify(cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang")
}

# GO enrichment for all DEGs (WD_cold vs SC_cold)
go_cold <- perform_go_enrichment(
  degs = df_degs_cold$ENTREZID,
  universe = data_wd_cold.vs.sc_cold$ENTREZID,
  orgdb = org.Mm.eg.db
)
write.xlsx(go_cold@result, "./reports/GO_DEGs_All_WD_cold_vs_SC_cold.xlsx")

# GO enrichment for upregulated DEGs
go_cold_up <- perform_go_enrichment(
  degs = df_degs_cold_up$ENTREZID,
  universe = data_wd_cold.vs.sc_cold$ENTREZID,
  orgdb = org.Mm.eg.db
)
write.xlsx(go_cold_up@result, "./reports/GO_DEGs_Up_WD_cold_vs_SC_cold.xlsx")

# GO enrichment for downregulated DEGs
go_cold_down <- perform_go_enrichment(
  degs = df_degs_cold_down$ENTREZID,
  universe = data_wd_cold.vs.sc_cold$ENTREZID,
  orgdb = org.Mm.eg.db
)
write.xlsx(go_cold_down@result, "./reports/GO_DEGs_Down_WD_cold_vs_SC_cold.xlsx")

# GO enrichment for interaction effect (Temperature x Diet)
go_temp_to_diet <- perform_go_enrichment(
  degs = df_degs_temp_to_diet$ENTREZID,
  universe = data_interaction_temp_to_diet$ENTREZID,
  orgdb = org.Mm.eg.db,
  ont = "all",
  p_cutoff = 0.1
)
write.xlsx(go_temp_to_diet@result, "./reports/GO_DEGs_All_Temperature_to_Diet.xlsx")

##########################################
# Step 2: KEGG Enrichment Analysis for Interaction Effect
##########################################

# KEGG enrichment analysis
kegg_temp_to_diet <- df_degs_temp_to_diet$ENTREZID %>%
  clusterProfiler::enrichKEGG(
    organism = "mmu",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
write.xlsx(kegg_temp_to_diet@result, "./reports/KEGG_DEGs_Temperature_to_Diet.xlsx")

# Convert KEGG ENTREZIDs to gene symbols
kegg_temp_to_diet@result <- convert_entrez_to_name(
  kegg_result = kegg_temp_to_diet@result,
  mapping_df = data_interaction_temp_to_diet,
  entrez_col = "geneID",
  id_col = "ENTREZID",
  name_col = "SYMBOL"
)

##########################################
# Step 3: GSEA for GO and KEGG Pathways
##########################################

# Prepare gene list for GSEA
gsea_genelist <- df_degs_temp_to_diet$log2FoldChange
names(gsea_genelist) <- df_degs_temp_to_diet$SYMBOL
gsea_genelist <- sort(gsea_genelist, decreasing = TRUE)

# Remove duplicates and NA values
gsea_genelist <- gsea_genelist[!is.na(names(gsea_genelist))]
gsea_genelist <- gsea_genelist[!duplicated(names(gsea_genelist))]

# Generate TERM2GENE for GO
gsea_gmt_go <- format_genes_to_term2gene(
  ID = go_temp_to_diet@result$ID,
  Description = go_temp_to_diet@result$Description,
  Genes = go_temp_to_diet@result$geneID,
  output_file = "./reports/GSEA/GO_temperature_to_diet.gmt"
)

# GSEA for GO pathways
gsea_go_temp_to_diet <- GSEA(
  geneList = gsea_genelist,
  TERM2GENE = gsea_gmt_go,
  pvalueCutoff = 0.6,
  minGSSize = 1,
  maxGSSize = 1000,
  eps = 1e-10,
  pAdjustMethod = "BH"
)
write.xlsx(gsea_go_temp_to_diet@result, "./reports/GSEA/GSEA_Result_GO_temperature_to_diet.xlsx")

# Generate TERM2GENE for KEGG
gsea_gmt_kegg <- format_genes_to_term2gene(
  ID = kegg_temp_to_diet@result$ID,
  Description = kegg_temp_to_diet@result$Description,
  Genes = kegg_temp_to_diet@result$geneID,
  output_file = "./reports/GSEA/KEGG_temperature_to_diet.gmt"
)

# GSEA for KEGG pathways
gsea_kegg_temp_to_diet <- GSEA(
  geneList = gsea_genelist,
  TERM2GENE = gsea_gmt_kegg,
  pvalueCutoff = 0.1,
  minGSSize = 3,
  maxGSSize = 1000,
  eps = 0,
  pAdjustMethod = "BH"
)
write.xlsx(gsea_kegg_temp_to_diet@result, "./reports/GSEA/GSEA_Result_KEGG_temperature_to_diet.xlsx")

##########################################
# Step 4: Combine and Save GSEA Results
##########################################

# Combine GO and KEGG GSEA results
gsea_combined_results <- rbind(
  as.data.frame(gsea_go_temp_to_diet@result),
  as.data.frame(gsea_kegg_temp_to_diet@result)
)
write.xlsx(gsea_combined_results, "./reports/GSEA/Combined_GSEA_Results.xlsx")

##########################################
# Step 5: Extract Gene Lists and Construct TF-Gene Network
##########################################

# Extract specific GO terms for TF analysis
tf_genes <- go_temp_to_diet@result %>%
  filter(ID == "GO_TERM_ID") %>% # Replace GO_TERM_ID with the desired term
  separate_rows(geneID, sep = "/") %>%
  select(geneID)

# Merge with DEG data to include log2FoldChange
tf_genes <- tf_genes %>%
  merge(data_interaction_temp_to_diet, by.x = "geneID", by.y = "ENTREZID") %>%
  select(SYMBOL, log2FoldChange)

# Save TF gene list
write.csv(tf_genes, "./reports/TRRUST/genelist.csv")

# Construct TF-Gene network
tf_enrichment <- read.xlsx("./reports/TRRUST/TRRUST_results_GO.xlsx", sheet = 1)

# Create edge list for network
edge_list <- tf_enrichment %>%
  separate_rows(`List.of.overlapped.genes`, sep = ",") %>%
  rename(TargetGene = `List.of.overlapped.genes`) %>%
  select(`Key.TF`, TargetGene)

# Create node attributes
tf_nodes <- tf_enrichment %>%
  select(`Key.TF`, `P.value`) %>%
  rename(name = `Key.TF`, attribute = `P.value`) %>%
  mutate(type = "TF")

gene_nodes <- tf_genes %>%
  rename(name = SYMBOL, attribute = log2FoldChange) %>%
  mutate(type = "Gene")

nodes <- bind_rows(tf_nodes, gene_nodes)

# Create network graph
network <- tbl_graph(nodes = nodes, edges = edge_list, directed = FALSE)

# Save network data
write.xlsx(as_data_frame(network, "edges"), "./reports/TRRUST/relationship.xlsx")
write.xlsx(as_data_frame(network, "nodes"), "./reports/TRRUST/TF_TargetGene_Network.xlsx")