library(Seurat)
library(ggplot2)
library(openxlsx)
library(snakecase)

root_dir <- "~/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/integrate_seurat_samples/")
dir.create(out_dir, showWarnings = F)

# read in metadata
metadata <- read.xlsx(paste0(root_dir, "HumanMeninges_Metadata_121125.xlsx"))

# get the samples
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


# add in metadata
metadata$Sample_name <- gsub("\\-", "", metadata$Sample_name)

all(metadata$Sample_name %in% seurat_merged$orig.ident)

metadata$condition <- to_snake_case(paste0(metadata$Group, "-", metadata$Age_Group))

metadata$Sample_name <- factor(as.character(metadata$Sample_name),
                               levels=metadata$Sample_name)

# get the sample metadata
sample_metadata <- seurat_merged@meta.data
sample_metadata$cell_id <- rownames(sample_metadata)

sample_metadata$Sample_name <- factor(sample_metadata$orig.ident,
                                      levels=levels(metadata$Sample_name))

sample_metadata <- merge(sample_metadata,
                         metadata,
                         by="Sample_name")
# correct the order
rownames(sample_metadata) <- sample_metadata$cell_id

sample_metadata <- sample_metadata[colnames(seurat_merged),]

# add to metadata in seurat object
seurat_merged@meta.data <- sample_metadata

# make some plots
DimPlot(seurat_merged, reduction="umap.harmony", 
        group.by= "harmony_clusters",
        split.by = "Sample_name",
        ncol=2, label = T)
ggsave(paste0(out_dir, "per_sample_umap.png"), width=7, height=9)

DimPlot(seurat_merged, reduction="umap.harmony", 
        group.by= "harmony_clusters",
        split.by = "condition",
        ncol=3, label = T)
ggsave(paste0(out_dir, "per_condition_umap.png"), width=9, height=4)

# save the integrated seurat
SaveSeuratRds(seurat_merged, paste0(out_dir, "all_samples.integrated_seurat.RDS"))

table(seurat_merged$orig.ident, seurat_merged$harmony_clusters)

# plot of the counts
cell_cluster_counts <- as.data.frame(table(sample_metadata$Sample_name,
                             sample_metadata$harmony_clusters))
cell_total_counts <- as.data.frame(table(sample_metadata$Sample_name))


# merge em
cell_cluster_counts <- merge(cell_cluster_counts,
                             cell_total_counts, by="Var1",
                             suffixes=c(".cluster",".total"))

cell_cluster_counts$frac <- (cell_cluster_counts$Freq.cluster / 
  cell_cluster_counts$Freq.total) * 100

cell_cluster_counts$Var2

ggplot(cell_cluster_counts,
       aes(x=Var2, y=frac, fill=Var1)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  labs(x="Harmony Cluster", y="Percentage of Cells in Cluster", fill=NULL)
ggsave(paste0(out_dir, "cluster_cell_fractions_per_sample.barplots.png"), width=10, height=5)


ggplot(cell_cluster_counts,
       aes(x=Var2, y=Freq.cluster, fill=Var1)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  labs(x="Harmony Cluster", y="Count of Cells in Cluster", fill=NULL)
ggsave(paste0(out_dir, "cluster_cell_counts_per_sample.barplots.png"), width=10, height=5)

ggplot(cell_total_counts,
       aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  guides(fill="none") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x=element_text())
ggsave(paste0(out_dir, "cell_counts_per_sample.barplots.png"), width=5, height=5)



