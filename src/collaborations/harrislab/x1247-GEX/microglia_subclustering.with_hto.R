library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(openxlsx)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/microglia_subclustering.with_hto/")
dir.create(out_dir, showWarnings = F)

# get seurat object
seu_obj <- readRDS(paste0(root_dir, "results/cell_type_naming.with_hto/seurat.cell_named.RDS"))

total_metadata <- seu_obj@meta.data

cluster2cell <- unique(total_metadata[,c("seurat_clusters","cell_name")])
cluster2cell_microglia <- cluster2cell[grepl("Microglia", cluster2cell$cell_name),]

microglia_cells <- rownames(total_metadata[total_metadata$seurat_clusters %in%
                                             cluster2cell_microglia$seurat_clusters,])

# subset
subset_seu <- subset(seu_obj, cells = microglia_cells)

# Doublet removal --------------------------------------------------------------

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Microglia Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")


max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "microglia_clusters")


# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

# do some doublet finding
# run doublet finder
subset_sce <- as.SingleCellExperiment(subset_seu)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
subset_sce <- scDblFinder(subset_sce, clusters="microglia_clusters", BPPARAM=bp)

# add doublet calls to seurat object

subset_seu$scDblFinder.class <- subset_sce$scDblFinder.class

VlnPlot(subset_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "microglia.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(subset_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "microglia.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(subset_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "microglia.doublet_output.RDS"))

# remove doublets, redo clustering
subset_seu <- subset(subset_seu, scDblFinder.class == "singlet")


# redo sc transform and harmony
subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Microglia Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "microglia_clusters", resolution = 0.8)

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.microglia_pca")

DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "microglia_clusters",
        label=T) 
ggsave(paste0(out_dir, "microglia_subcluster.umap.png"), width=5, height=4)

DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= "microglia_clusters",
        label=T, split.by = "condition", ncol=2) 
ggsave(paste0(out_dir, "microglia_sublcuster.condition.umap.png"), width=6, height=5)

# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

Idents(subset_seu) <- "microglia_clusters"

all_markers <- FindAllMarkers(subset_seu)

# save results
SaveSeuratRds(subset_seu, paste0(out_dir, "subset_microglia.seurat.RDS"))
saveRDS(all_markers, paste0(out_dir, "all_markers.microglia.RDS"))
saveRDS(subset_seu@meta.data, paste0(out_dir, "subset_microglia.metadata.RDS"))




