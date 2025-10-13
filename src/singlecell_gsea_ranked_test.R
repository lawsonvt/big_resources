# Example of how to do ranked GSEA on single cell data
# Two approaches:
# Using Seurat FindMarkers() for simple comparisons
# Using pseudobulk + DESeq2 for comparisons involving replicates

library(Seurat)
library(SeuratData)
library(fgsea)  # Fast GSEA implementation
library(msigdbr)  # For gene sets
library(dplyr)
library(ggplot2)
library(DESeq2)
library(cowplot)

# get data for doing analysis
InstallData("ifnb")

LoadData("ifnb")
# Update old Seurat object to accommodate new features
ifnb <- UpdateSeuratObject(ifnb)


# get gene sets (Hallmark gene set)
# Get gene sets from MSigDB (e.g., Hallmark pathways)
gene_sets <- msigdbr(species = "Homo sapiens", collection = "H")

# Convert to list format for fgsea
pathways <- gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Or use other collections:
# "C2" for curated gene sets
# "C5" for GO terms
# "C8" for cell type signatures

# first way to do it, use FindMarkers on a per cell basis
Idents(ifnb) <- "seurat_annotations"

# focus on one cell
ifnb_cd14 <- subset(ifnb, idents = "CD14 Mono")

# label cells with conditions
Idents(ifnb_cd14) <- "stim"
table(Idents(ifnb_cd14))

# find differentially expressed genes between markers
markers <- FindMarkers(
  ifnb_cd14,
  ident.1 = "STIM",
  ident.2 = "CTRL",
  test.use = "wilcox",  # or "MAST", "DESeq2"
  logfc.threshold = 0,  # Get all genes
  min.pct = 0.1
)
markers$gene_name <- rownames(markers)

# find minimum non-zero p-value
min_nz_pval <- min(markers$p_val[markers$p_val >0])

# create "stat" value for ranking
markers$stat <- -log10(pmax(markers$p_val, min_nz_pval)) * sign(markers$avg_log2FC)

# Create ranked gene list (by log fold change or stat)
# Important: remove NAs and sort
ranked_genes <- markers %>%
  filter(!is.na(stat) &
           !is.infinite(stat)) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = gene_name)

# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500
)

fgsea_results %>%
  arrange(padj) %>%
  head(20)


# Table plot of top pathways
topPathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  head(20) %>%
  pull(pathway)

plotGseaTable(
  pathways[topPathways],
  ranked_genes,
  fgsea_results,
  gseaParam = 0.5
)


# Plot enrichment for top pathway
p1 <- plotEnrichment(
  pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], 
  ranked_genes
) + labs(title = "IFN-γ Response")

p2 <- plotEnrichment(
  pathways[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], 
  ranked_genes
) + labs(title = "Oxidative Phosphorylation")

plot_grid(p1, p2)


# Pseudobulk analysis, featuring DESeq2

# aggregate cells to create "pseudobulk" 
aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), 
                                      return.seurat = TRUE)
Idents(aggregate_ifnb) <- "seurat_annotations"
agg_ifnb_cdmono <- subset(aggregate_ifnb, idents=c("CD14 Mono",
                                                   "CD16 Mono"))

# create metadata for samples
sample_metadata <- agg_ifnb_cdmono@meta.data %>%
  select(stim, seurat_annotations, orig.ident) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$orig.ident

# create DESeq data set
dds <- DESeqDataSetFromMatrix(
  countData = agg_ifnb_cdmono[["RNA"]]$counts,
  colData = sample_metadata,
  design = ~ stim
)

# run DESeq
dds <- DESeq(dds)

# pull out results
res <- results(dds, contrast = c("stim", "STIM", "CTRL"))

# reformat for extracting
res <- as.data.frame(res)
res$gene_name <- rownames(res)

# Create ranked list
ranked_genes <- res %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = gene_name)

# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500
)

# Table plot of top pathways
topPathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  pull(pathway)

plotGseaTable(
  pathways[topPathways],
  ranked_genes,
  fgsea_results,
  gseaParam = 0.5
)

# Plot enrichment for top pathway
plotEnrichment(
  pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], 
  ranked_genes
) + labs(title = "IFN-γ Response")

# compare top pathway to bottom pathway
p1 <- plotEnrichment(
  pathways[["HALLMARK_INFLAMMATORY_RESPONSE"]], 
  ranked_genes
) + labs(title = "Inflammatory Response") + ylim(-0.2, 0.6)

p2 <- plotEnrichment(
  pathways[["HALLMARK_XENOBIOTIC_METABOLISM"]], 
  ranked_genes
) + labs(title = "Xenobiotic Metabolism") + ylim(-0.2, 0.6)

plot_grid(p1, p2)
