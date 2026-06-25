library(Seurat)
library(openxlsx)
library(ggplot2)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/find_all_markers/")
dir.create(out_dir, showWarnings = F, recursive = T)


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

all_markers <- FindAllMarkers(int_seu)

saveRDS(all_markers, file=paste0(out_dir, "all_markers.RDS"))




