library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(openxlsx)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)
library(dplyr)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/subclustering_tcell/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

# group the conditions
metadata <- seu_obj@meta.data

metadata$sample <- factor(metadata$orig.ident,
                          levels=c("X1161_1NEG",
                                   "X1176-WT1",
                                   "X1177-WT2",
                                   "X1162_2Het",
                                   "X1178-KO1",
                                   "X1179-KO2"))

metadata$condition <- "WT"

metadata[metadata$sample %in% c("X1162_2Het",
                                "X1178-KO1",
                                "X1179-KO2"),]$condition <- "KO"
metadata$condition <- factor(metadata$condition,
                             levels=c("WT","KO"))

# add to seurat object
seu_obj$sample <- metadata$sample
seu_obj$condition <- metadata$condition

# Filtering main set down to T cells -------------------------------------------

# load in immune subclustering metadata
immune_metadata <- readRDS(paste0(root_dir, "results/immune_cell_split/immune_cell_subclustering/subset_immune.metadata.RDS"))


# get cluster numbers
immune_cluster_raw <- read.xlsx(paste0(root_dir, "immune cell classification.xlsx"))

# reformat
immune_cluster_celltypes <- lapply(immune_cluster_raw$Cell.type, function(cell) {
  
  clusters <- unlist(strsplit(immune_cluster_raw[immune_cluster_raw$Cell.type == cell,]$Cluster, ","))
  
  return(data.frame(celltype=cell,
                    cluster=clusters))
  
})
immune_cluster_celltypes <- bind_rows(immune_cluster_celltypes)

# filter dwon to tcells
tcell_clusters <- immune_cluster_celltypes[immune_cluster_celltypes$celltype == "T cell",]

tcell_metadata <- immune_metadata[immune_metadata$immune_clusters %in% tcell_clusters$cluster,]

# filter down and do subclustring
subset_seu <- subset(seu_obj, cells=rownames(tcell_metadata))

# Doublet removal --------------------------------------------------------------

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="T Cell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "tcell_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")


max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "tcell_clusters")


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
subset_sce <- scDblFinder(subset_sce, clusters="tcell_clusters", BPPARAM=bp)

# add doublet calls to seurat object

subset_seu$scDblFinder.class <- subset_sce$scDblFinder.class

VlnPlot(subset_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "tcell.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(subset_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "tcell.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(subset_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "tcell.doublet_output.RDS"))

# remove doublets, redo clustering
subset_seu <- subset(subset_seu, scDblFinder.class == "singlet")

# Re-cluster post doublet removal

# redo sc transform and harmony
subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="sample")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="tcell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "tcell_clusters")

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T) 
ggsave(paste0(out_dir, "tcell_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T, split.by="condition") 
ggsave(paste0(out_dir, "tcell_subcluster.per_condition.umap.png"), width=10, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, split.by="sample", ncol=3) 
ggsave(paste0(out_dir, "tcell_subcluster.per_sample.umap.png"), width=14, height=7)

# clusters 2 and 6 are exclusive to one samples, remove and redo

# Remove single samples clusters and redo clustering ---------------------------

subset_seu <- subset(subset_seu, subset = !tcell_clusters %in% c("4","11", "12"))

table(subset_seu@meta.data$sample,
      subset_seu@meta.data$tcell_clusters)


# redo sc transform and harmony
subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="sample")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="tcell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu,  cluster.name = "tcell_clusters")

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T) 
ggsave(paste0(out_dir, "tcell_subcluster.post_filter.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T, split.by="condition") 
ggsave(paste0(out_dir, "tcell_subcluster.per_condition.post_filter.umap.png"), width=10, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, split.by="sample", ncol=3) 
ggsave(paste0(out_dir, "tcell_subcluster.per_sample.post_filter.umap.png"), width=14, height=7)

# determine markers and save results
# find markers!
# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

Idents(subset_seu) <- "tcell_clusters"

all_markers <- FindAllMarkers(subset_seu)

# save results
SaveSeuratRds(subset_seu, paste0(out_dir, "subset_tcell.seurat.RDS"))
saveRDS(all_markers, paste0(out_dir, "all_markers.tcell.RDS"))
saveRDS(subset_seu@meta.data, paste0(out_dir, "subset_tcell.metadata.RDS"))

# bar plots
tcell_meta <- subset_seu@meta.data

ggplot(tcell_meta,
       aes(x=sample)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=4,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "tcell_cluster_counts.sample_bar_plot.png"), width=10, height=8)

ggplot(tcell_meta,
       aes(x=condition)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=4,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "tcell_cluster_counts.condition_bar_plot.png"), width=10, height=8)









