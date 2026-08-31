library(Seurat)


root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/find_all_markers.immune_cells/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/immune_cell_subclustering/subset_immune.seurat.RDS"))

Idents(seu_obj) <- "immune_clusters"

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

all_markers <- FindAllMarkers(seu_obj)

saveRDS(all_markers, file=paste0(out_dir, "all_markers.RDS"))


