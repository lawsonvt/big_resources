library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)
library(cowplot)
library(gtools)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/immune_cell.pca_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in nonimmune seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/immune_cell_subclustering/subset_immune.seurat.RDS"))

metadata <- seu_obj@meta.data

# add in cluster word to clusters
levels(metadata$immune_clusters) <- paste0("cluster", levels(metadata$immune_clusters))

seu_obj$immune_clusters <- metadata$immune_clusters

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

# make total pseudo clusters

total_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("sample")
)


# Extract the count matrix
total_pseudo_matrix <- total_pseudo$RNA
# fix sample names
colnames(total_pseudo_matrix) <- gsub("\\-1NEG", "_1NEG", colnames(total_pseudo_matrix))
colnames(total_pseudo_matrix) <- gsub("\\-2Het", "_2Het", colnames(total_pseudo_matrix))

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

total_pseudo_matrix <- total_pseudo_matrix[,rownames(sample_metadata)]

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = total_pseudo_matrix,
  colData = sample_metadata,
  design = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind=F)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA, all cells") +
  geom_text_repel(aes(label=sample)) +
  theme_bw()
ggsave(paste0(out_dir, "total.pca_plot.png"), width=6, height=5)

# do the same analysis on a per cluster basis

cluster_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("immune_clusters","sample")
)

# Extract the count matrix
cluster_pseudo_matrix <- cluster_pseudo$RNA
# fix sample names
colnames(cluster_pseudo_matrix) <- gsub("\\-1NEG", "_1NEG", colnames(cluster_pseudo_matrix))
colnames(cluster_pseudo_matrix) <- gsub("\\-2Het", "_2Het", colnames(cluster_pseudo_matrix))


clusters <- levels(metadata$immune_clusters)

# do just the first 20 clusters, since later ones have missing cell counts
clusters <- clusters[1:20]

cluster_dir <- paste0(out_dir, "cluster_pcas/")
dir.create(cluster_dir, showWarnings = F)

pca_plot_list <- lapply(clusters, function(cluster) {
  
  print(cluster)
  # get counts
  counts <- cluster_pseudo_matrix[,grepl(paste0(cluster, "_"), colnames(cluster_pseudo_matrix))]
  
  # fix column names
  colnames(counts) <- gsub(paste0(cluster, "_"), "", colnames(counts))
  
  all(colnames(counts) %in% rownames(sample_metadata)) &
    all(rownames(sample_metadata) %in% colnames(counts))
  
  #counts <- counts[,rownames(sample_metadata)]
  
  counts_meta <- sample_metadata[colnames(counts),]
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = counts_meta,
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  
  vsd <- vst(dds, blind=F)
  
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste0("PCA, ", cluster)) +
    geom_text_repel(aes(label=sample), size=3, max.overlaps = 50) +
    theme_bw() +
    theme(legend.position = "none")
  p
  ggsave(paste0(cluster_dir, cluster, ".pca_plot.png"), width=6, height=5)
    
  return(p)
})

plot_grid(plotlist = pca_plot_list, ncol=5)
ggsave(paste0(out_dir, "all_clusters_pca_plot.png"), width=14, height=12)

# look at the PCAs after having removed NEG/HET

samples <- levels(seu_obj$sample)
samples <- grep("(NEG)|(Het)", samples, value=T, invert=T)

subset_seu <- subset(seu_obj, subset = sample %in% samples)

total_pseudo <- AggregateExpression(
  subset_seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("sample")
)


# Extract the count matrix
total_pseudo_matrix <- total_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- subset_seu@meta.data %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

total_pseudo_matrix <- total_pseudo_matrix[,rownames(sample_metadata)]

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = total_pseudo_matrix,
  colData = sample_metadata,
  design = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind=F)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA, all cells") +
  geom_text_repel(aes(label=sample)) +
  theme_bw()
ggsave(paste0(out_dir, "total.filtered.pca_plot.png"), width=6, height=5)


# do the same analysis on a per cluster basis

cluster_pseudo <- AggregateExpression(
  subset_seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("immune_clusters","sample")
)

# Extract the count matrix
cluster_pseudo_matrix <- cluster_pseudo$RNA


clusters <- levels(metadata$immune_clusters)

# do just the first 20 clusters, since later ones have missing cell counts
clusters <- clusters[1:20]

cluster_dir <- paste0(out_dir, "cluster_pcas_filtered/")
dir.create(cluster_dir, showWarnings = F)

pca_plot_list <- lapply(clusters, function(cluster) {
  
  print(cluster)
  # get counts
  counts <- cluster_pseudo_matrix[,grepl(paste0(cluster, "_"), colnames(cluster_pseudo_matrix))]
  
  # fix column names
  colnames(counts) <- gsub(paste0(cluster, "_"), "", colnames(counts))
  
  all(colnames(counts) %in% rownames(sample_metadata)) &
    all(rownames(sample_metadata) %in% colnames(counts))
  
  #counts <- counts[,rownames(sample_metadata)]
  
  counts_meta <- sample_metadata[colnames(counts),]
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = counts_meta,
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  
  vsd <- vst(dds, blind=F)
  
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste0("PCA, ", cluster)) +
    geom_text_repel(aes(label=sample), size=3, max.overlaps = 50) +
    theme_bw() +
    theme(legend.position = "none")
  p
  ggsave(paste0(cluster_dir, cluster, ".pca_plot.png"), width=6, height=5)
  
  return(p)
})

plot_grid(plotlist = pca_plot_list, ncol=5)
ggsave(paste0(out_dir, "all_clusters_filtered_pca_plot.png"), width=18, height=10)


