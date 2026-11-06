library(Seurat)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(dplyr)

root_dir <- "~/Documents/projects/gaultierlab/shreyas_gujar/snRNASeq/"

# test comment

out_dir <- paste0(root_dir, "results/seurat_qc_samples/")

dir.create(out_dir, showWarnings = F, recursive = T)

# read in all the samples
h5_files <- list("O1"=paste0(root_dir, "cellranger_outputs/O1/filtered_feature_bc_matrix.h5"),
                 "O2"=paste0(root_dir, "cellranger_outputs/O2/filtered_feature_bc_matrix.h5"),
                 "O3"=paste0(root_dir, "cellranger_outputs/O3/filtered_feature_bc_matrix.h5"),
                 "S1"=paste0(root_dir, "cellranger_outputs/S1/filtered_feature_bc_matrix.h5"),
                 "S2"=paste0(root_dir, "cellranger_outputs/S2/filtered_feature_bc_matrix.h5"),
                 "S3"=paste0(root_dir, "cellranger_outputs/S3/filtered_feature_bc_matrix.h5"))

print("Read in samples ...")

# read in samples
samples <- lapply(names(h5_files), function(sample_name) {
  
  print(sample_name)
  
  sample <- Read10X_h5(h5_files[[sample_name]])
  
  return(CreateSeuratObject(counts=sample, project=sample_name, min.cells = 3, min.features = 200))
  
})
names(samples) <- names(h5_files)

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

filter_counts_list <- list()

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  
  # make a data frame of 
  filter_counts <- data.frame(type=c("Raw Cell Count",
                                     paste0("nFeature_RNA > ", feature_min),
                                     paste0("nFeature_RNA < ", feature_max),
                                     paste0("percent.mt < ", percent_mt_max)),
                              count=c(nrow(sample@meta.data),
                                      sum(sample@meta.data$nFeature_RNA > feature_min),
                                      sum(sample@meta.data$nFeature_RNA < feature_max),
                                      sum(sample@meta.data$percent.mt < percent_mt_max)))
  
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
  
  # filter bar counts
  filter_counts <- rbind(filter_counts,
                         data.frame(type="Filtered Count",
                                    count=nrow(sample@meta.data)))
  
  filter_counts$type <- factor(filter_counts$type,
                               levels=filter_counts$type)
  
  ggplot(filter_counts,
         aes(x=type,
             y=count)) +
    geom_bar(stat="identity",
             color="black",
             fill="maroon") +
    geom_text(aes(label=count),
              vjust=-0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=35, hjust=1)) +
    labs(x=NULL, y="Cell Count")
  ggsave(paste0(plot_dir, sample_name, ".filter_parameter_counts.png"), width=4, height=5)
  
  # save filter counts for a combined plot
  
  filter_counts$sample <- sample_name
  
  filter_counts_list <<- c(filter_counts_list,
                          list(filter_counts))
  
  return(sample)
})

filter_counts_df <- bind_rows(filter_counts_list)

ggplot(filter_counts_df,
       aes(x=type,
           y=count)) +
  geom_bar(stat="identity",
           color="black",
           fill="maroon") +
  geom_text(aes(label=count),
            vjust=-0.3) +
  facet_wrap(~ sample, ncol=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(plot_dir, "total.filter_parameter_counts.png"), width=9, height=7)

# sctransform and PCA elbow plot

plot_dir <- paste0(out_dir, "3_init_pca_plots/")
dir.create(plot_dir, showWarnings = F)

print("Initial PCA plots ...")

samples <- lapply(samples, function(sample) {
  
  sample_name <- sample@project.name
  
  print(sample_name)
  
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

# doublet finding!
# code adapted from tutorial 
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

# redo pca elbow plots, to ensure PCA number is still good

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

# clustering

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

