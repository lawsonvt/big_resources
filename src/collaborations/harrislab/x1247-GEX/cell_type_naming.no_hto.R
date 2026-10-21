library(Seurat)
library(SeuratObject)
library(ggplot2)
library(openxlsx)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/cell_type_naming.no_hto/")
dir.create(out_dir, showWarnings = F)

seu_obj <- readRDS(paste0(root_dir, "results/seurat_qc_samples.no_hto/seurat.post_qc.RDS"))

cell_annotations <- read.xlsx(paste0(root_dir, "results/cassia_annotate.no_hto/celltype_labels.xlsx"))

metadata <- seu_obj@meta.data

metadata$cell_id <- rownames(metadata)

# merge in cell names
metadata <- merge(metadata,
                  cell_annotations,
                  by.x="seurat_clusters",
                  by.y="Cluster.ID")

# fix order
rownames(metadata) <- metadata$cell_id

metadata <- metadata[colnames(seu_obj),]

seu_obj@meta.data <- metadata

DimPlot(seu_obj, reduction="umap", group.by=c("seurat_clusters", "cell_name"),
        label=T)
ggsave(paste0(out_dir, "clusters.umap.png"), width=18, height=8)

SaveSeuratRds(seu_obj, paste0(out_dir, "seurat.cell_named.RDS"))


