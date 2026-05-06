library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(openxlsx)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/subclustering_tcell/")
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

# subset down to tcell
subset_seu <- subset(seu_obj, subset = cell_type %in% c("T Cell",
                                                        "Cytotoxic Treg-like Cells",
                                                        "Treg"))

# clean up and free memory
rm(seu_obj)
gc()

# initial umap cluster
DimPlot(subset_seu, reduction="umap.harmony", group.by= "cell_type",
        label=F, split.by = "sample", ncol=3) 
ggsave(paste0(out_dir, "init_tcell_clusters.samples.umap.png"), width=10, height=7)

DimPlot(subset_seu, reduction="umap.harmony", group.by= "cell_type",
        label=F, split.by = "condition") 
ggsave(paste0(out_dir, "init_tcell_clusters.condition.umap.png"), width=8, height=5)

# lets do some doublet detection!

# determine the optimal number of dimensions through an elbow plot

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="T Cell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "tcell_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

# 20 dimensions makes sense for this dataset as well

max_pc_dim <- 20

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


# there is going to need to be some filtering, so the first set of plots are going to
# be put into a pre_filter folder

pre_dir <- paste0(out_dir, "pre_filter/")
dir.create(pre_dir, showWarnings = F)

# determine the optimal number of dimensions through an elbow plot

# redo sc transform and harmony
subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="sample")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="tcell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(pre_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "tcell_clusters")

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T) 
ggsave(paste0(pre_dir, "tcell_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T, split.by="condition") 
ggsave(paste0(pre_dir, "tcell_subcluster.per_condition.umap.png"), width=10, height=5)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, split.by="sample", ncol=3) 
ggsave(paste0(pre_dir, "tcell_subcluster.per_sample.umap.png"), width=14, height=7)

DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "cell_type",
        label=F, split.by="sample", ncol=3) 
ggsave(paste0(pre_dir, "tcell_subcluster.init_cluster.per_sample.umap.png"), width=14, height=7)

# DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= "cell_type",
#         label=F) 
# ggsave(paste0(out_dir, "tcell_orig_cell_type.umap.png"), width=8, height=6)

# cluster counts
tcell_meta <- subset_seu@meta.data

ggplot(tcell_meta,
       aes(x=sample)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(pre_dir, "tcell_cluster_counts.sample_bar_plot.png"), width=10, height=8)

ggplot(tcell_meta,
       aes(x=condition)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(pre_dir, "tcell_cluster_counts.condition_bar_plot.png"), width=10, height=8)

# try filtering the clump that only belongs to one sample

sample_clump <- c("3","13","18")

# filter down
subset_subset_seu <- subset(subset_seu, subset = !tcell_clusters %in% sample_clump)

# redo sc transform and harmony
subset_subset_seu <- SCTransform(subset_subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_subset_seu <- RunPCA(subset_subset_seu, npcs = 50)

subset_subset_seu <- RunHarmony(subset_subset_seu, group.by.vars="sample")

# inspect elbow plot
ElbowPlot(subset_subset_seu, ndims=50) + 
  labs(title="tcell Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
subset_subset_seu <- FindNeighbors(subset_subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_subset_seu <- FindClusters(subset_subset_seu, cluster.name = "tcell_clusters")

# create umap
subset_subset_seu <- RunUMAP(subset_subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(subset_subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T) 
ggsave(paste0(out_dir, "tcell_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_subset_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",
        label=T, split.by="condition") 
ggsave(paste0(out_dir, "tcell_subcluster.per_condition.umap.png"), width=10, height=5)

DimPlot(subset_subset_seu, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, split.by="sample", ncol=3) 
ggsave(paste0(out_dir, "tcell_subcluster.per_sample.umap.png"), width=14, height=7)

DimPlot(subset_subset_seu, reduction="umap.tcell_pca", group.by= "cell_type",
        label=F, split.by="sample", ncol=3) 
ggsave(paste0(out_dir, "tcell_subcluster.init_cluster.per_sample.umap.png"), width=14, height=7)

# cluster counts
tcell_meta <- subset_subset_seu@meta.data

ggplot(tcell_meta,
       aes(x=sample)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "tcell_cluster_counts.sample_bar_plot.png"), width=10, height=8)

ggplot(tcell_meta,
       aes(x=condition)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ tcell_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "tcell_cluster_counts.condition_bar_plot.png"), width=10, height=8)


# find markers!

clusters <- sort(unique(subset_subset_seu$tcell_clusters))

# set assay to RNA and join layers
DefaultAssay(subset_subset_seu) <- "RNA"
subset_subset_seu <- JoinLayers(subset_subset_seu)

subset_subset_seu <- NormalizeData(subset_subset_seu)
subset_subset_seu <- FindVariableFeatures(subset_subset_seu)
subset_subset_seu <- ScaleData(subset_subset_seu)

cluster_markers <- lapply(clusters, function(cluster) {
  
  markers <- FindMarkers(subset_subset_seu,
                         ident.1 = cluster,
                         test.use = "wilcox")
  markers <- rownames_to_column(markers, var = "gene")
  
  markers$pct_delta <- markers$pct.1 - markers$pct.2
  
  markers$cluster <- cluster
  
  return(markers[order(markers$pct_delta, decreasing = T),])
  
})

# save to file
write.xlsx(cluster_markers,
           paste0(out_dir, "cluster_markers.xlsx"),
           colWidths="auto")


# output top 10 for each cluster
cluster_dir <- paste0(out_dir, "cluster_marker_dot_plots/")

dir.create(cluster_dir, showWarnings = F)

Idents(subset_subset_seu) <- "tcell_clusters"

for (data in cluster_markers) {
  
  top_markers <- data$gene[1:10]
  
  cluster_name <- unique(data$cluster)
  
  DotPlot(subset_subset_seu, features = top_markers) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=cluster_name) 
  ggsave(paste0(cluster_dir, "cluster", cluster_name, ".top_markers_dots.png"),
         width=6, height=10, bg="white")
  
}

SaveSeuratRds(subset_subset_seu, paste0(out_dir, "subset_tcell.seurat.RDS"))
