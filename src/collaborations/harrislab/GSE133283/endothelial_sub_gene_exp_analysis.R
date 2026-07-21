library(openxlsx)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(plyr)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"


out_dir <- paste0(root_dir, "results/endothelial_sub_gene_exp_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in Seurat object

endo_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/endothelial_subclustering/endothelial_subcluster.seurat.RDS"))

# ensure correct clusters
Idents(endo_seu) <- "endothelial_clusters"

# exam the genes of interest
group1_genes <- c("Il1r1",
                  "Il1rap",
                  "Tnfrsf1a", 
                  "Ifngr1",
                  "Ifngr2")

group1_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group1_genes)
ggsave(paste0(out_dir, "group1_violin_plots.png"), width=9, height=5)

DotPlot(endo_seu, features=group1_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group1_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group1_genes, reduction="umap.endothelial_pca", ncol=3)
ggsave(paste0(out_dir, "group1_umap.png"), width=12, height=8)

FeaturePlot(endo_seu, features = group1_genes, 
            reduction = "umap.endothelial_pca",
            split.by="condition")
ggsave(paste0(out_dir, "group1_umap.per_condition.png"), width=12, height=14)

# look at the next set of genes of interest
group2_genes <- c("B2m",
                  "Tap1",
                  "Cxcl10",
                  "H2-Ab1")


group2_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group2_genes, ncol=2)
ggsave(paste0(out_dir, "group2_violin_plots.png"),  width=8, height=7)

DotPlot(endo_seu, features=group2_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group2_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group2_genes, reduction="umap.endothelial_pca", ncol=2)
ggsave(paste0(out_dir, "group2_umap.png"), width=8, height=7)

FeaturePlot(endo_seu, features = group2_genes, 
            reduction = "umap.endothelial_pca",
            split.by="condition")
ggsave(paste0(out_dir, "group2_umap.per_condition.png"), width=12, height=12)

# group 3!
group3_genes <- c("Stat1","Ciita", "Icam1", "Vcam1")

group3_genes %in% rownames(endo_seu)

VlnPlot(endo_seu, features = group3_genes, ncol=2)
ggsave(paste0(out_dir, "group3_violin_plots.png"),  width=8, height=7)

DotPlot(endo_seu, features=group3_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "group3_dot_plots.png"), width=6, height=5, bg="white")


FeaturePlot(endo_seu, features = group3_genes, reduction="umap.endothelial_pca", ncol=2)
ggsave(paste0(out_dir, "group3_umap.png"), width=8, height=7)

FeaturePlot(endo_seu, features = group3_genes, 
            reduction = "umap.endothelial_pca",
            split.by="condition")
ggsave(paste0(out_dir, "group3_umap.per_condition.png"), width=12, height=12)


# determine the markers that separate the clusters

endo_markers <- readRDS(paste0(root_dir,
                               "results/endothelial_subclustering/all_markers.RDS"))
endo_markers$pct_diff <- endo_markers$pct.1 - endo_markers$pct.2
endo_markers <- endo_markers[order(endo_markers$cluster),]

endo_markers_top <- endo_markers[endo_markers$pct_diff > 0.4,]

DotPlot(endo_seu, features=unique(endo_markers_top$gene)) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "top_pct_genes.dot_plot.png"), width=14, height=6, bg="white")

endo_clusters <- sort(unique(endo_markers$cluster))

endo_markers <- endo_markers[order(endo_markers$pct_diff, decreasing = T),]

#pull out top endomarkers
top_endo_markers <- lapply(endo_clusters, function(x) {
  
  head(endo_markers[endo_markers$cluster == x &
                      endo_markers$pct_diff > 0,])
  
})
top_endo_markers <- bind_rows(top_endo_markers)


DotPlot(endo_seu, features=unique(top_endo_markers$gene)) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "top_rank_genes.dot_plot.png"), width=16, height=7, bg="white")

# filter down to adult and EAE conditions

endo_seu_ae <- subset(endo_seu, subset = condition %in% c("Adult", "EAE"))

VlnPlot(endo_seu_ae, features = group1_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group1_violin_plots.condition_split.png"), width=6, height=11)


VlnPlot(endo_seu_ae, features = group2_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group2_violin_plots.condition_split.png"), width=6, height=10)

VlnPlot(endo_seu_ae, features = group3_genes, split.by = "condition", ncol=1)
ggsave(paste0(out_dir, "group3_violin_plots.condition_split.png"), width=6, height=10)


