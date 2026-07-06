library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/immune_cell_subclustering/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir, "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# list of immune clusters
subset_clusters <- c("5","8","11","16")

# subset down to immune cells
subset_seu <- subset(seu_obj, subset = harmony_clusters %in% subset_clusters)

# clean up and free memory
rm(seu_obj)
gc()

subset_seu <- RunPCA(subset_seu, npcs = 50)

# calculated this wrong in initial QC ...
subset_seu[["percent.mt"]] <- PercentageFeatureSet(subset_seu, pattern = "^MT-")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Immune Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "immune_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "immune_clusters")

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
subset_sce <- scDblFinder(subset_sce, clusters="immune_clusters", BPPARAM=SerialParam(RNGseed = 1234),
                          samples = "Sample_name")

# add doublet calls to seurat object

subset_seu$scDblFinder.class <- subset_sce$scDblFinder.class

VlnPlot(subset_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "immune.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(subset_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "immune.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(subset_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "immune.doublet_output.RDS"))

# remove doublets, redo clustering
subset_seu <- subset(subset_seu, scDblFinder.class == "singlet")

# set Assay back to SCT for clustering
#DefaultAssay(subset_seu) <- "SCT"

subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="Sample_name", theta = 4)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Immune Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "immune_subset.post_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

# 20 dimensions makes sense for this dataset as well

max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "harmony")
subset_seu <- FindClusters(subset_seu, cluster.name = "immune_clusters", resolution = 0.8)

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="harmony", reduction.name="umap.immune_pca")

DimPlot(subset_seu, reduction="umap.immune_pca", group.by= "immune_clusters",
        label=T) 
ggsave(paste0(out_dir, "immune_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.immune_pca", group.by= "immune_clusters",
        label=T, split.by = "condition", ncol = 3) 
ggsave(paste0(out_dir, "immune_subcluster.per_sample.umap.png"), width=12, height=5)

# cluster counts
subset_meta <- subset_seu@meta.data

ggplot(subset_meta,
       aes(x=Sample_name)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ immune_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "immune_cluster_counts.sample_bar_plot.png"), width=10, height=8)

ggplot(subset_meta,
       aes(x=condition)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ immune_clusters, ncol=6,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "immune_cluster_counts.condition_bar_plot.png"), width=10, height=8)

clusters <- sort(unique(subset_seu$immune_clusters))

# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

cluster_markers <- lapply(clusters, function(cluster) {
  
  markers <- FindMarkers(subset_seu,
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

Idents(subset_seu) <- "immune_clusters"

for (data in cluster_markers) {
  
  top_markers <- data$gene[1:10]
  
  cluster_name <- unique(data$cluster)
  
  DotPlot(subset_seu, features = top_markers) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=cluster_name) 
  ggsave(paste0(cluster_dir, "cluster", cluster_name, ".top_markers_dots.png"),
         width=6, height=10, bg="white")
  
}

saveRDS(subset_seu@meta.data, paste0(out_dir, "subset_immune.metadata.RDS"))

SaveSeuratRds(subset_seu, paste0(out_dir, "subset_immune.seurat.RDS"))

all_markers <- FindAllMarkers(subset_seu)

saveRDS(all_markers, paste0(out_dir, "all_markers.RDS"))

