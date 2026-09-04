library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/subclustering_opc.round2/")
dir.create(out_dir, showWarnings = F)

# load seurat data
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/nonimmune_cell_subclustering/subset_nonimmune.seurat.RDS"))

metadata <- seu_obj@meta.data

# subset down to probable OPC
subset_seu <- subset(seu_obj, subset = nonimmune_clusters %in% c("2","10","18","19","20"))

# plot to ensure correct
DimPlot(subset_seu, reduction="umap.nonimmune_pca", 
        group.by= "nonimmune_clusters",
        label=T)
# there are some cells that seem very split from the UMAP, maybe investigate later?

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="OPC Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "opc_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 15

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "opc_clusters")

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
subset_sce <- scDblFinder(subset_sce, clusters="opc_clusters", BPPARAM=bp)

# add doublet calls to seurat object

subset_seu$scDblFinder.class <- subset_sce$scDblFinder.class

VlnPlot(subset_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "opc.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(subset_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "opc.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(subset_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "opc.doublet_output.RDS"))

# remove doublets, redo clustering
subset_seu <- subset(subset_seu, scDblFinder.class == "singlet")

# set Assay back to SCT for clustering
#DefaultAssay(subset_seu) <- "SCT"

subset_seu <- SCTransform(subset_seu, vars.to.regress = c("percent.mt"), verbose = F)

subset_seu <- RunPCA(subset_seu, npcs = 50)

subset_seu <- RunHarmony(subset_seu, group.by.vars="sample")

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="OPC Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "opc_subset.post_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

# 20 dimensions makes sense for this dataset 

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "opc_clusters", resolution = 0.1)


# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.opc_pca")

DimPlot(subset_seu, reduction="umap.opc_pca", group.by= "opc_clusters",
        label=T) 
ggsave(paste0(out_dir, "opc_subcluster.umap.png"), width=7, height=5)

DimPlot(subset_seu, reduction="umap.opc_pca", group.by= "opc_clusters",
        label=T, split.by = "condition") 
ggsave(paste0(out_dir, "opc_subcluster.per_condition.umap.png"), width=10, height=5)

# find markers!
# set assay to RNA and join layers
DefaultAssay(subset_seu) <- "RNA"
subset_seu <- JoinLayers(subset_seu)

subset_seu <- NormalizeData(subset_seu)
subset_seu <- FindVariableFeatures(subset_seu)
subset_seu <- ScaleData(subset_seu)

Idents(subset_seu) <- "opc_clusters"

all_markers <- FindAllMarkers(subset_seu)

# save results
SaveSeuratRds(subset_seu, paste0(out_dir, "subset_opc.seurat.RDS"))
saveRDS(all_markers, paste0(out_dir, "all_markers.opc.RDS"))

# cluster counts
subset_meta <- subset_seu@meta.data

ggplot(subset_meta,
       aes(x=sample)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ opc_clusters, ncol=2,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "opc_cluster_counts.sample_bar_plot.png"), width=7, height=5)

ggplot(subset_meta,
       aes(x=condition)) +
  geom_bar(color="black", fill="grey") +
  facet_wrap(~ opc_clusters, ncol=2,
             scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "opc_cluster_counts.condition_bar_plot.png"), width=7, height=5)

# export markers
markers_list <- lapply(unique(all_markers$cluster), function(cluster) {
  
  data <- all_markers[all_markers$cluster == cluster,]
  # reorder columns
  data <- data[,c("gene", setdiff(colnames(data), "gene"))]
  data$delta_pct <- data$pct.1 - data$pct.2
  
  return(data)
})
names(markers_list) <- unique(all_markers$cluster)

write.xlsx(markers_list, file=paste0(out_dir, "cluster_markers.xlsx"), colWidths="auto")





