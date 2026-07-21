library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tibble)
library(scDblFinder)
library(BiocParallel)
library(openxlsx)
library(harmony)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/endothelial_subclustering/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in Seurat object

seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# total metadata
total_meta <- seu_obj@meta.data

total_meta$sample_id <- total_meta$orig.ident
total_meta$condition <- gsub("\\-[0-9]+", "", total_meta$sample_id)
total_meta$condition <- factor(total_meta$condition)

# add in cluster word to clusters
levels(total_meta$harmony_clusters) <- paste0("cluster", levels(total_meta$harmony_clusters))

seu_obj$harmony_clusters <- total_meta$harmony_clusters

seu_obj$sample_id <- total_meta$sample_id
seu_obj$condition <- total_meta$condition

# plots!
ggplot(total_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")

ggplot(total_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ harmony_clusters, scale="free", ncol=6) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")

# identify endothelial clusters and subset
endothelial_clusters <- c("cluster0",
                          "cluster2",
                          "cluster4",
                          "cluster6",
                          "cluster11",
                          "cluster16",
                          "cluster18")


endo_seu <- subset(seu_obj, subset = harmony_clusters %in% endothelial_clusters)

endo_meta <- endo_seu@meta.data


ggplot(endo_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "endothelial_subcluster.sample_cell_counts.png"),
       width=7, height=5)

ggplot(endo_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ harmony_clusters, scale="free", ncol=3) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "endothelial_total_cluster.per_cluster_sample_counts.png"),
       width=8, height=7)

# try reclustering

endo_seu <- RunPCA(endo_seu, npcs = 50)


# inspect elbow plot
ElbowPlot(endo_seu, ndims=50) + 
  labs(title="Endothelial Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "endothelial_subset.pre_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

max_pc_dim <- 20

# cluster the harmonized data
endo_seu <- FindNeighbors(endo_seu, dims = 1:max_pc_dim, reduction = "pca")
endo_seu <- FindClusters(endo_seu, cluster.name = "endothelial_clusters")

# set assay to RNA and join layers
DefaultAssay(endo_seu) <- "RNA"
endo_seu <- JoinLayers(endo_seu)

endo_seu <- NormalizeData(endo_seu)
endo_seu <- FindVariableFeatures(endo_seu)
endo_seu <- ScaleData(endo_seu)

# do some doublet finding
# run doublet finder
endo_sce <- as.SingleCellExperiment(endo_seu)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
endo_sce <- scDblFinder(endo_sce, clusters="endothelial_clusters", BPPARAM=SerialParam(RNGseed = 1234),
                          samples = "sample_id")

endo_seu$scDblFinder.class <- endo_sce$scDblFinder.class

VlnPlot(endo_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "endothelial.doublet_qc_plots.png"), width=9, height=6)


doublet_counts <- as.data.frame(table(endo_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "endothelial.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(endo_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "endothelial.doublet_output.RDS"))

# remove doublets, redo clustering
endo_seu <- subset(endo_seu, scDblFinder.class == "singlet")

# set Assay back to SCT for clustering
#DefaultAssay(endo_seu) <- "SCT"

endo_seu <- SCTransform(endo_seu, vars.to.regress = c("percent.mt"), verbose = F)

endo_seu <- RunPCA(endo_seu, npcs = 50)

endo_seu <- RunHarmony(endo_seu, group.by.vars="sample_id", theta = 2)

# inspect elbow plot
ElbowPlot(endo_seu, ndims=50) + 
  labs(title="Endothelial Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "endothelial_subset.post_doublet.pca_elbow_plot.png"), width=6, height=5, bg="white")

# 20 dimensions makes sense for this dataset as well

max_pc_dim <- 20

# cluster the harmonized data
endo_seu <- FindNeighbors(endo_seu, dims = 1:max_pc_dim, reduction = "harmony")
endo_seu <- FindClusters(endo_seu, cluster.name = "endothelial_clusters", resolution = 0.8)

# create umap
endo_seu <- RunUMAP(endo_seu, dims = 1:max_pc_dim, reduction="harmony", reduction.name="umap.endothelial_pca")


DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "endothelial_clusters",
        label=T, label.box = T) 
ggsave(paste0(out_dir, "endothelial_subcluster.umap.png"), width=7, height=5)

DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "endothelial_clusters",
        label=T, split.by = "condition", ncol = 2) 
ggsave(paste0(out_dir, "endothelial_subcluster.per_sample.umap.png"), width=8, height=7)

# make more plots
endo_meta <- endo_seu@meta.data


ggplot(endo_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ endothelial_clusters, scale="free", ncol=5) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "endothelial_subcluster.per_cluster_sample_counts.png"),
       width=12, height=8)

ggplot(endo_meta,
       aes(x=endothelial_clusters, 
           fill=condition)) +
  geom_bar() +
  theme_bw() +
  labs(x="Endothelial Subclusters", y="Cell Count")
ggsave(paste0(out_dir, "endothelial_subcluster.cluster_counts.png"),
       width=8, height=6)

# filter out low count clusters / non diverse clusters

cluster_condition_counts <- as.data.frame(table(endo_meta$endothelial_clusters,
                                                endo_meta$condition))
colnames(cluster_condition_counts) <- c("endothelial_clusters",
                                        "condition",
                                        "cell_count")

min_cell_count <- 100

cluster_condition_counts_f <- cluster_condition_counts[cluster_condition_counts$cell_count >= min_cell_count,]

keep_clusters <- sort(unique(cluster_condition_counts_f$endothelial_clusters[duplicated(cluster_condition_counts_f$endothelial_clusters)]))

endo_seu <- subset(endo_seu, subset = endothelial_clusters %in% keep_clusters)



DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "endothelial_clusters",
        label=T, label.box = T) 
ggsave(paste0(out_dir, "filtered_endothelial_subcluster.umap.png"), width=7, height=5)

DimPlot(endo_seu, reduction="umap.endothelial_pca", group.by= "endothelial_clusters",
        label=T, split.by = "condition", ncol = 2) 
ggsave(paste0(out_dir, "filtered_endothelial_subcluster.per_sample.umap.png"), width=8, height=7)

# make more plots
endo_meta <- endo_seu@meta.data

ggplot(endo_meta,
       aes(x=sample_id,
           fill=condition)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ endothelial_clusters, scale="free", ncol=4) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "filtered_endothelial_subcluster.per_cluster_sample_counts.png"),
       width=11, height=7)

ggplot(endo_meta,
       aes(x=endothelial_clusters, 
           fill=condition)) +
  geom_bar() +
  theme_bw() +
  labs(x="Endothelial Subclusters", y="Cell Count")
ggsave(paste0(out_dir, "filtered_endothelial_subcluster.cluster_counts.png"),
       width=7, height=5)

# find all markers

Idents(endo_seu) <- "endothelial_clusters"

# set assay to RNA and join layers
DefaultAssay(endo_seu) <- "RNA"
endo_seu <- JoinLayers(endo_seu)

endo_seu <- NormalizeData(endo_seu)
endo_seu <- FindVariableFeatures(endo_seu)
endo_seu <- ScaleData(endo_seu)

all_markers <- FindAllMarkers(endo_seu)

saveRDS(all_markers, file=paste0(out_dir, "all_markers.RDS"))

# save file
SaveSeuratRds(endo_seu, paste0(out_dir, "endothelial_subcluster.seurat.RDS"))













