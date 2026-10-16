library(Seurat)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)

root_dir <- "~/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/seurat_qc_samples/")


h5_file <- paste0(root_dir, "cellranger_outputs/filtered_feature_bc_matrix.h5")