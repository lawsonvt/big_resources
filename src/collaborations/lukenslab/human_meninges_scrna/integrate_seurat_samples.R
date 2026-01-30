library(Seurat)
library(ggplot2)
library(scRepertoire)

root_dir <- "~/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/integrate_seurat_samples/")
dir.create(out_dir, showWarnings = F)


sample_dir <- paste0(root_dir, "results/seurat_qc_samples/")
# read in seurat objects

seurat_files <- list.files(path=sample_dir, 
                           pattern=".qc_processed.seurat.RDS",
                           full.names = T)

# load in files
seurat_list <- lapply(seurat_files, function(file) {
  
  LoadSeuratRds(file)
  
})

# pull in names
names(seurat_list) <- unlist(lapply(seurat_list, function(x) {x@project.name}))

# merge them ...
seurat_merged <- merge(x=seurat_list[[1]],
                       y=seurat_list[-1],
                       add.cell.ids=names(seurat_list),
                       merge.data=T) # merge normalized data as well

# save the merged seurat
SaveSeuratRds(seurat_merged, paste0(out_dir, "all_samples.merged_seurat.RDS"))


# quiet the TCR genes from being variable features
# seurat_merged <- quietTCRgenes(seurat_merged)

# no normalization required, since we merge the normalized data as well
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)

# inspect elbow plot
ElbowPlot(seurat_merged, ndims=30) + 
  labs(title="Merged Samples (not integrated)") +
  scale_x_continuous(breaks=seq(0,30,5)) +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
ggsave(paste0(out_dir, "merge_seurat.pca_elbow_plot.png"), width=6, height=5, bg="white")


# after inspecting the elbowplot
max_pc_dim <- 20

# first do some clustering with the merged data
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:max_pc_dim, reduction = "pca")
seurat_merged <- FindClusters(seurat_merged, cluster.name = "unintegrated_clusters")

# plot out a umap which will showcase the batch effect
seurat_merged <- RunUMAP(seurat_merged, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.unintegrated")

DimPlot(seurat_merged, reduction="umap.unintegrated", group.by= "orig.ident")
ggsave(paste0(out_dir, "merge_seurat.umap.png"), width=8, height=6)

# INTEGRATE
seurat_merged <- IntegrateLayers(seurat_merged,
                                 method=HarmonyIntegration,
                                 orig.reduction = "pca",
                                 new.reduction = "harmony",
                                 verbose = F)

# cluster the harmonized data
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:max_pc_dim, reduction = "harmony")
seurat_merged <- FindClusters(seurat_merged, cluster.name = "harmony_clusters")

# create umap
seurat_merged <- RunUMAP(seurat_merged, dims = 1:max_pc_dim, reduction="harmony", reduction.name="umap.harmony")

# plot it
DimPlot(seurat_merged, reduction="umap.harmony", group.by= c("orig.ident","harmony_clusters"))
ggsave(paste0(out_dir, "integrated_seurat.umap.png"), width=14, height=6)

# save the integrated seurat
SaveSeuratRds(seurat_merged, paste0(out_dir, "all_samples.integrated_seurat.RDS"))

table(seurat_merged$orig.ident, seurat_merged$harmony_clusters)





