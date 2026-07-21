library(openxlsx)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(plyr)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"


out_dir <- paste0(root_dir, "results/endothelial_sub_contamination_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in Seurat object

endo_seu <- LoadSeuratRds(paste0(root_dir,
                                 "results/endothelial_subclustering/endothelial_subcluster.seurat.RDS"))

# ensure correct clusters
Idents(endo_seu) <- "endothelial_clusters"

# suggested by Claude
markers <- c("Cldn5","Pecam1","Flt1", # endothelial
             "Aqp4","Gfap", # astrocyte
             "Rbfox3","Snap25", # neuron
             "Cx3cr1","P2ry12", # microglia
             "Col1a", "Lum") # fibroblast

DotPlot(endo_seu, features = markers)






