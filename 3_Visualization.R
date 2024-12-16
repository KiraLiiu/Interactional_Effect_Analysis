##########################################
# RNA-seq Acylcarnitine Project - Visualization
# Author: Kira Liu
# Date: 2024-11-18
# Purpose: Create visualizations for DEG analysis, including volcano plots, GO/KEGG bar plots,
#          cnetplots, GSEA plots, and transcription factor (TF) network visualizations.
##########################################

#### Load Required Libraries ####
library(VennDiagram)      # Venn diagrams
library(ComplexHeatmap)   # Heatmap visualizations
library(circlize)         # Color settings for circular plots
library(RColorBrewer)     # Color palettes
library(ggplot2)          # Data visualization
library(ggrepel)          # Labeling significant genes
library(viridis)          # Advanced color maps
library(stringr)          # String manipulation
library(enrichplot)       # Enrichment analysis visualizations
library(tidyverse)        # Data manipulation and visualization
library(tidygraph)        # Network analysis
library(ggraph)           # Graph plotting

##########################################
# Volcano Plots for DEGs
##########################################

# Add Significance column for WD_cold vs SC_cold
data_wd_cold.vs.sc_cold <- data_wd_cold.vs.sc_cold %>%
  mutate(Significance = case_when(
    log2FoldChange > log2fc_threshold & padj < pvalue_threshold ~ "Upregulated",
    log2FoldChange < -log2fc_threshold & padj < pvalue_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Volcano plot for WD_cold vs SC_cold
ggplot(data_wd_cold.vs.sc_cold, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  scale_color_manual(name = "DEGs",
                     values = c("Upregulated" = "#de2d26", "Downregulated" = "#08519c", "Not Significant" = "lightgrey")) +
  labs(x = "log2FC", y = "-log10(p.adjust)", title = "DEGs: WD_cold vs SC_cold") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1)
  ) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 45)) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "solid") +
  geom_text_repel(
    data = filter(data_wd_cold.vs.sc_cold, Significance != "Not Significant"),
    aes(label = SYMBOL),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = 'grey50'
  )

# Save volcano plot
ggsave("./plots/volcanoplot_DEGs_WD_cold.vs.SC_cold.png", width = 8, height = 6, dpi = 300)

##########################################
# GO Enrichment Bar Plot
##########################################

# Combine GO enrichment results
go_cold_up@result$Category <- "Up"
go_cold_down@result$Category <- "Down"
go_cold_combined <- rbind(go_cold_up@result, go_cold_down@result)

# Select top 5 terms by adjusted p-value for each category
top5 <- go_cold_combined %>%
  group_by(Category) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

# Reverse order of descriptions for plotting
top5$Description <- factor(top5$Description, levels = rev(top5$Description))

# Create GO bar plot
p <- ggplot(top5, aes(x = Count, y = Description, fill = Category)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0.5)) +
  scale_fill_manual(values = alpha(c('#fa9a6a', '#6d79c9'), 0.66)) +
  labs(title = "Top GO Terms for WD_cold vs SC_cold", x = "Gene Count", y = "GO Term") +
  theme(
    plot.title = element_text(size = 22, face = 'bold', hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13)
  )

# Save GO bar plot
ggsave("./plots/GO_Terms_Barplot_WD_cold.vs.SC_cold.png", p, width = 9, height = 7, dpi = 300)

##########################################
# Circular Network Plot for GO Terms
##########################################

# Generate gene list for cnetplot
enrich_geneList <- data_interaction_temp_to_diet$log2FoldChange
names(enrich_geneList) <- data_interaction_temp_to_diet$SYMBOL

# Create circular network plot
cnetplot(go_temp_to_diet, foldChange = enrich_geneList, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Circular Network Plot: GO Terms")

# Save circular plot
ggsave("./plots/GO_Circular_Network_Plot.png", width = 8, height = 6, dpi = 300)

##########################################
# GSEA Visualization
##########################################

# GSEA plot for GO terms
gseaplot2(
  gsea_go_temp_to_diet,
  geneSetID = c("GO:0001216"),
  title = "GSEA for GO Function: WD_cold vs SC_cold",
  color = "#394B96",
  pvalue_table = TRUE,
  base_size = 11
)

# Save GSEA plot
ggsave("./plots/GSEA/GSEA_GO_WD_cold.vs.SC_cold.png", width = 16, height = 8, dpi = 300)

##########################################
# TF-Gene Network Visualization
##########################################

# Combine nodes and edges
nodes <- bind_rows(tf_nodes, gene_nodes)
network <- tbl_graph(nodes = nodes, edges = edge_list, directed = FALSE)

# Plot network
ggraph(network, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5), color = "grey") +
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_graph() +
  ggtitle("TF-Gene Regulatory Network")

# Save TF-Gene network plot
ggsave("./plots/TF_Gene_Network.png", width = 10, height = 8, dpi = 300)