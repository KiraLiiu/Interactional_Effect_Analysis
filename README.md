# RNA-seq Analysis Pipeline

## Overview

This repository contains an RNA-seq differential expression analysis (DEA) pipeline implemented in R. The analysis focuses on studying the main effects (diet and temperature) and interaction effects (temperature’s influence on dietary responses) in a mouse model. The pipeline includes data preprocessing, DEA, enrichment analysis (GO and KEGG), GSEA, and advanced visualization (e.g., volcano plots, bar plots, and network plots).

## Analysis Workflow

1. Data Preprocessing
	•	Raw count data and metadata are loaded from an Excel file.
	•	Preprocessing steps:
	•	Ensure count data are integers.
	•	Create combined condition labels (e.g., WD_cold, SC_cold).
	•	Filter out low-expressed genes (row sum > 10).

2. Differential Expression Analysis (DEA)
	•	Conducted using DESeq2:
	•	Compare gene expression across conditions for main effects:
	•	Dietary effect: WD_cold vs SC_cold, WD_RT vs SC_RT
	•	Temperature effect: WD_cold vs WD_RT, SC_cold vs SC_RT
	•	Evaluate the interaction effect: Temperature x Diet.
	•	Differentially expressed genes (DEGs) are identified based on:
	•	|log2FoldChange| > threshold
	•	Adjusted p-value (padj) < 0.05

3. Functional Enrichment Analysis
	•	GO Enrichment (Biological Processes, BP):
	•	GO terms enriched in all DEGs, upregulated DEGs, and downregulated DEGs are identified using clusterProfiler.
	•	Simplified results using semantic similarity with a cutoff (0.5).
	•	KEGG Pathway Enrichment:
	•	Enrichment analysis performed using enrichKEGG.
	•	Pathway terms are mapped to gene symbols for interpretability.

4. Gene Set Enrichment Analysis (GSEA)
	•	GSEA is conducted for both GO terms and KEGG pathways:
	•	Input: Ranked gene list (log2FoldChange).
	•	Gene sets formatted into TERM2GENE for GSEA.
	•	Results include enrichment scores (ES) and visualized plots.

5. Visualization

5.1 Volcano Plots
	•	Visualize DEGs by plotting:
	•	log2FoldChange (x-axis)
	•	-log10(p-adjusted value) (y-axis)
	•	Annotate significant DEGs.

5.2 GO/KEGG Bar Plots
	•	Display top enriched terms based on adjusted p-values.
	•	Group terms into upregulated and downregulated categories.

5.3 Circular Network Plots (Cnetplots)
	•	Visualize relationships between enriched GO terms and associated genes.
	•	Circular layout and color-coded edges for better readability.

5.4 GSEA Plots
	•	Display GSEA results, including enrichment scores and p-value tables.

5.5 TF-Gene Regulatory Network
	•	Construct and visualize transcription factor (TF) networks:
	•	Nodes: Transcription factors (TFs) and target genes.
	•	Edges: Regulatory relationships.
	•	Output: Network graphs for exploration.

## Key Features
	•	Interaction Analysis: Captures interaction effects between diet and temperature using DESeq2.
	•	Comprehensive Functional Insights:
	•	GO/KEGG enrichment provides biological context for DEGs.
	•	GSEA highlights broader pathway-level perturbations.
	•	Advanced Visualizations:
	•	Publication-ready plots for volcano plots, bar plots, and networks.
	•	Reproducibility:
	•	Modular code organized by workflow stages.
	•	Results saved in well-defined directories.
