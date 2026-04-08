library(Seurat)
library(SeuratObject)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

# output directory
out_dir <- paste0(root_dir, "results/split_and_integrate/")
dir.create(out_dir, showWarnings = F)

# read in seurat
seu_obj <- LoadSeuratRds(paste0(root_dir, "data/Data.Combined_clustering.rds"))
metadata <- seu_obj@meta.data

# since QC was already done, maybe just split and integrate?

seu_list <- SplitObject(seu_obj, split.by = "orig.ident")

# slim down the Seurat objects
seu_list <- lapply(seu_list, function(x) {
  
  DietSeurat(x, layers="counts")
  
})

# rename the layers
seu_list[[1]]@assays$RNA@layers$counts <- seu_list[[1]]@assays$RNA@layers$counts.1
seu_list[[1]]@assays$RNA@layers$counts.1 <- NULL

seu_list[[2]]@assays$RNA@layers$counts <- seu_list[[2]]@assays$RNA@layers$counts.2
seu_list[[2]]@assays$RNA@layers$counts.2 <- NULL

seu_list[[3]]@assays$RNA@layers$counts <- seu_list[[3]]@assays$RNA@layers$counts.3
seu_list[[3]]@assays$RNA@layers$counts.3 <- NULL

seu_list[[4]]@assays$RNA@layers$counts <- seu_list[[4]]@assays$RNA@layers$counts.4
seu_list[[4]]@assays$RNA@layers$counts.4 <- NULL

# SC transform
seu_list <- lapply(seu_list, function(x) {
  SCTransform(x, vars.to.regress = "percent.mt")
})

# lets integrate!
seu_merged <- merge(x=seu_list[[1]],
                       y=seu_list[-1],
                       add.cell.ids=NULL, # they are already there
                       merge.data=T) # merge normalized data as well

# save the merged seurat
SaveSeuratRds(seu_merged, paste0(out_dir, "all_samples.merged_seurat.RDS"))


# no normalization required, since we merge the normalized data as well
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- ScaleData(seu_merged)
seu_merged <- RunPCA(seu_merged)

ElbowPlot(seu_merged, ndims=30) + 
  labs(title="Merged Samples (not integrated)") +
  scale_x_continuous(breaks=seq(0,30,5)) +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
ggsave(paste0(out_dir, "merge_seurat.pca_elbow_plot.png"), width=6, height=5, bg="white")

# after inspecting the elbowplot
max_pc_dim <- 20

# first do some clustering with the merged data
seu_merged <- FindNeighbors(seu_merged, dims = 1:max_pc_dim, reduction = "pca")
seu_merged <- FindClusters(seu_merged, cluster.name = "unintegrated_clusters")

# plot out a umap which will showcase the batch effect
seu_merged <- RunUMAP(seu_merged, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.unintegrated")

# unintegrated umap
DimPlot(seu_merged, reduction="umap.unintegrated", group.by= "orig.ident")
ggsave(paste0(out_dir, "merge_seurat.umap.png"), width=8, height=6)

# INTEGRATE
seu_merged <- IntegrateLayers(seu_merged,
                                 method=HarmonyIntegration,
                                 orig.reduction = "pca",
                                 new.reduction = "harmony",
                                 verbose = F)

# cluster the harmonized data
seu_merged <- FindNeighbors(seu_merged, dims = 1:max_pc_dim, reduction = "harmony")
seu_merged <- FindClusters(seu_merged, cluster.name = "harmony_clusters")

# create umap
seu_merged <- RunUMAP(seu_merged, dims = 1:max_pc_dim, reduction="harmony", reduction.name="umap.harmony")

# plot it
DimPlot(seu_merged, reduction="umap.harmony", group.by= c("orig.ident","harmony_clusters"))
ggsave(paste0(out_dir, "integrated_seurat.umap.png"), width=14, height=6)

DimPlot(seu_merged, reduction="umap.harmony", group.by= c("MULTI_ID","harmony_clusters"))

DimPlot(seu_merged, reduction="umap.harmony", group.by="harmony_clusters",
        split.by="MULTI_ID", ncol=2)
ggsave(paste0(out_dir, "integrated_seurat.umap_per_sample.png"), width=8, height=10)


# save the integrated file
SaveSeuratRds(seu_merged, paste0(out_dir, "all_samples.integrated_seurat.RDS"))


