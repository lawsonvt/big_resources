library(Seurat)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

# output directory
out_dir <- paste0(root_dir, "results/explore_seurat/")

dir.create(out_dir, recursive = T, showWarnings = F)

seu <- LoadSeuratRds(paste0(root_dir, "data/Data.Combined_clustering.rds"))

metadata <- seu@meta.data

# can i join the data?
seu <- JoinLayers(seu)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)


# make some dim plots with the UMAP
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label=T) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "cluster_umap.png"), width=7, height=5)

DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label=T,
        split.by = "MULTI_ID", ncol=2) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "cluster_umap.sample_split.png"), width=8, height=10)


DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label=T,
        split.by = "orig.ident", ncol=2) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "cluster_umap.library_split.png"), width=8, height=7)

# add in cell type information
metadata$cell_type <- as.character(metadata$seurat_clusters)

metadata[metadata$seurat_clusters %in% c("5","6","7","13"),]$cell_type <- "Astrocytes"
metadata[metadata$seurat_clusters %in% c("4"),]$cell_type <- "Microglia"

seu$cell_type <- metadata$cell_type


DimPlot(seu, reduction = "umap", group.by = "cell_type", label=T) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "celltypes_umap.png"), width=7, height=5)


DimPlot(seu, reduction = "umap", group.by = "cell_type", label=T,
        split.by = "MULTI_ID", ncol=2) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "celltypes_umap.sample_split.png"), width=8, height=10)

DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label=T,
        split.by = "orig.ident", ncol=2) +
  labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "celltypes_umap.library_split.png"), width=8, height=7)