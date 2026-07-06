library(Seurat)
library(SeuratObject)
library(ggplot2)
library(openxlsx)

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/find_all_markers/")
dir.create(out_dir, showWarnings = F)

seu_obj <- LoadSeuratRds(paste0(root_dir, "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

Idents(seu_obj) <- "harmony_clusters"

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

all_markers <- FindAllMarkers(seu_obj)

saveRDS(all_markers, file=paste0(out_dir, "all_markers.RDS"))
