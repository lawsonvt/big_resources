library(Seurat)
library(ggplot2)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/Harris-AM-6617/"


out_dir <- paste0(root_dir, "results/marker_dot_plots/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

Idents(int_seu) <- "harmony_clusters"

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))

all_markers$pct_diff <- all_markers$pct.1 - all_markers$pct.2

sum(all_markers$pct_diff > 0.75)

table(all_markers[all_markers$pct_diff > 0.75,]$cluster)

all_markers <- all_markers[order(all_markers$pct_diff, decreasing = T),]

#pull out top allmarkers
top_all_markers <- lapply(sort(unique(all_markers$cluster)), function(x) {
  
  head(all_markers[all_markers$cluster == x &
                      all_markers$pct_diff > 0,], n=3)
  
})
top_all_markers <- bind_rows(top_all_markers)

DotPlot(int_seu, features=unique(top_all_markers$gene)) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "top_rank_genes.dot_plot.png"), width=16, height=7, bg="white")

# the point was to enrich for astrocyte markers, can we find any?

# source https://pmc.ncbi.nlm.nih.gov/articles/PMC9265979/#sec4-cells-11-02021
# source https://www.nature.com/articles/s41467-021-20892-3/figures/4

astro_markers <- c("Gfap",
                   "S100b",
                   "Aldh1l1",
                   "Slc1a2",
                   "Slc1a3",
                   "Sox9",
                   "Aqp4")

DotPlot(int_seu, features=astro_markers) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "astro_markers.dot_plot.png"), width=8, height=7, bg="white")




