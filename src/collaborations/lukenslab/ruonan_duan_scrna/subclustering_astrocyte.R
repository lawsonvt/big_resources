library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(openxlsx)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/subclustering_astrocyte/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

# subset down to astrocyte
subset_seu <- subset(seu_obj, subset = cell_type == "Astrocyte")

# clean up and free memory
rm(seu_obj)
gc()

# determine the optimal number of dimensions through an elbow plot

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Astrocyte Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "astrocyte_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 10

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "astrocyte_clusters")

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.astrocyte_pca")

DimPlot(subset_seu, reduction="umap.astrocyte_pca", group.by= "astrocyte_clusters",
        label=T) 
ggsave(paste0(out_dir, "astrocyte_subcluster.umap.png"), width=7, height=5)

# find markers!

clusters <- sort(unique(subset_seu$astrocyte_clusters))

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

Idents(subset_seu) <- "astrocyte_clusters"

for (data in cluster_markers) {
  
  top_markers <- data$gene[1:10]
  
  cluster_name <- unique(data$cluster)
  
  DotPlot(subset_seu, features = top_markers) + 
    RotatedAxis() + 
    labs(x=NULL, y=NULL, title=cluster_name) 
  ggsave(paste0(cluster_dir, "cluster", cluster_name, ".top_markers_dots.png"),
         width=6, height=10, bg="white")
  
}

SaveSeuratRds(subset_seu, paste0(out_dir, "subset_astrocyte.seurat.RDS"))
