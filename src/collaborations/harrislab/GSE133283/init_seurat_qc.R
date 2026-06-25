library(Seurat)
library(SeuratObject)
library(ggplot2)
library(readr)
library(reshape2)
library(scDblFinder)
library(BiocParallel)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/init_seurat_qc/")
dir.create(out_dir, showWarnings = F, recursive = T)

# get the files
count_files <- list.files(paste0(root_dir, "GSE133283_RAW/"),
                          "counts.tsv.gz",
                          full.names=T)

# get sample IDs
names(count_files) <- sapply(count_files, function(x) {
  unlist(strsplit(basename(x), "_"))[2]
})

# create seurat objects

samples <- lapply(names(count_files), function(sample_id) {
  
  count_file <- count_files[sample_id]
  
  counts <- read_tsv(count_file)
  
  counts_mat <- acast(counts, gene ~ cell, value.var = "count")
  
  counts_mat[is.na(counts_mat)] <- 0
  
  return(CreateSeuratObject(counts=counts_mat, project=sample_id, min.cells = 3, min.features = 200))
  
})
names(samples) <- names(count_files)

plot_dir <- paste0(out_dir, "1_pre_filter_qc_plots/")
dir.create(plot_dir, showWarnings = F)

print("Pre-filtering QC plots ...")

samples <- lapply(samples, function(sample) {
  sample_name <- sample@project.name
  
  print(sample_name)
  
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
  
  VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.4) 
  ggsave(paste0(plot_dir, sample_name, ".pre_filter_count_plot.png"), width=8, height=6)
  
  plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  ggsave(paste0(plot_dir, sample_name, ".pre_filter_feature_scatter.png"), width=12, height=6)
  
  return(sample)
})

# filter and replot
feature_min <- 200
feature_max <- 6000 # neuronal cells express MANY genes, so we can make this high
percent_mt_max <- 5

plot_dir <- paste0(out_dir, "2_post_filter_qc_plots/")
dir.create(plot_dir, showWarnings = F)

print("Filter and replot QC plots ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  
  # filter
  sample <- subset(sample, 
                   subset=nFeature_RNA > feature_min &
                     nFeature_RNA < feature_max &
                     percent.mt < 5)
  
  VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.4) 
  ggsave(paste0(plot_dir, sample_name, ".post_filter_count_plot.png"), width=8, height=6)
  
  plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  ggsave(paste0(plot_dir, sample_name, ".post_filter_feature_scatter.png"), width=12, height=6)
  
  return(sample)
})

# sctransform and PCA elbow plot

plot_dir <- paste0(out_dir, "3_init_pca_plots/")
dir.create(plot_dir, showWarnings = F)

print("Initial PCA plots ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  # quiet the TCR genes from being variable features
  #sample <- quietTCRgenes(sample)
  
  # regress out mitochondrial percentage for SC Transform
  sample <- SCTransform(sample, vars.to.regress = "percent.mt")
  
  sample <- RunPCA(sample)
  
  ElbowPlot(sample, ndims=30) + 
    labs(title=sample_name) +
    scale_x_continuous(breaks=seq(0,30,5)) +
    scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
  ggsave(paste0(plot_dir, sample_name, ".pca_elbow_plot.png"), width=6, height=5, bg="white")
  
  return(sample)
})

# it looks like the elbow based number for dimensions should be 20, at least for now

max_pc_dim <- 20

# doublet finding!
# code adated from tutorial 
# https://plger.github.io/scDblFinder/articles/scDblFinder.html

plot_dir <- paste0(out_dir, "4_doublet_detection/")
dir.create(plot_dir, showWarnings = F)

print("Doublet Detection and Removal ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  
  # find clusters, purely for doublet finding
  sample <- FindNeighbors(sample, dims = 1:max_pc_dim)
  sample <- FindClusters(sample)
  
  sce <- as.SingleCellExperiment(sample)
  
  bp <- MulticoreParam(3, RNGseed=1234) # equivalent to set seed, for reproducibility
  sce <- scDblFinder(sce, clusters="seurat_clusters", BPPARAM=bp)
  
  sample$scDblFinder.class <- sce$scDblFinder.class
  
  VlnPlot(sample, group.by="orig.ident", split.by = "scDblFinder.class",
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3, pt.size = 0) + theme(legend.position = 'right')
  ggsave(paste0(plot_dir, sample_name, ".doublet_qc_plots.png"), width=9, height=6)
  
  doublet_counts <- as.data.frame(table(sample$scDblFinder.class))
  
  ggplot(doublet_counts, 
         aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", fill="grey", color="black") +
    geom_text(aes(label=Freq), vjust=-0.4) +
    theme_bw() +
    labs(x=NULL, y="Cell Count", title=sample_name)
  ggsave(paste0(plot_dir, sample_name, ".doublet_counts.png"), width=4, height=5)
  
  dbl_meta <- as.data.frame(colData(sce))
  
  saveRDS(dbl_meta, file=paste0(plot_dir, sample_name, ".doublet_output.RDS"))
  
  # final step, remove the doublets
  
  return(subset(sample, scDblFinder.class == "singlet"))
  
})

plot_dir <- paste0(out_dir, "5_post_doublet_pca_plots/")
dir.create(plot_dir, showWarnings = F)

print("Post Doublet PCA plots ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  # quiet the TCR genes from being variable features
  #sample <- quietTCRgenes(sample)
  
  # regress out mitochondrial percentage for SC Transform
  sample <- SCTransform(sample, vars.to.regress = "percent.mt")
  
  sample <- RunPCA(sample)
  
  ElbowPlot(sample, ndims=30) + 
    labs(title=sample_name) +
    scale_x_continuous(breaks=seq(0,30,5)) +
    scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
  ggsave(paste0(plot_dir, sample_name, ".pca_elbow_plot.png"), width=6, height=5, bg="white")
  
  return(sample)
})

max_pc_dim <- 20

plot_dir <- paste0(out_dir, "6_clustering_umap/")
dir.create(plot_dir, showWarnings = F)

print("Clustering and UMAP ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  
  sample <- FindNeighbors(sample, dims=1:max_pc_dim)
  sample <- FindClusters(sample)
  
  sample <- RunUMAP(sample, dims=1:max_pc_dim)
  
  DimPlot(sample, reduction="umap")
  ggsave(paste0(plot_dir, sample_name, ".cluster_umap.png"), width=8, height=6)
  
  return(sample)
})

# output final samples

for (sample in samples) {
  
  sample_name <- sample@project.name
  
  SaveSeuratRds(sample, file=paste0(out_dir, sample_name, ".qc_processed.seurat.RDS"))
  
}

