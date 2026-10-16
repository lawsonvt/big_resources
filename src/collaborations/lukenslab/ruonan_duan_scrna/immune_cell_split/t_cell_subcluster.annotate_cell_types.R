library(Seurat)
library(SeuratObject)
library(openxlsx)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/t_cell_subcluster.annotate_cell_type/")
dir.create(out_dir, showWarnings = F)

# read in sublcutered T cell
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/immune_cell_split/subclustering_tcell/subset_tcell.seurat.RDS"))

# read in cell names
cluster_xref <- read.xlsx(paste0(root_dir, "immune_cell_split.tcell_subcluster.cell_types.xlsx"))

metadata <- seu_obj@meta.data

metadata$cell_id <- rownames(metadata)
metadata$cell_type <- NULL

cluster_xref$cluster <- factor(as.character(cluster_xref$cluster),
                               levels=levels(metadata$tcell_clusters))

# merge in t cell names
metadata <- merge(metadata,
                  cluster_xref,
                  by.x="tcell_clusters",
                  by.y="cluster")
rownames(metadata) <- metadata$cell_id

# fix order
metadata <- metadata[colnames(seu_obj),]

# add to seurat object
seu_obj@meta.data <- metadata

DimPlot(seu_obj, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, label.box=T) 
ggsave(paste0(out_dir, "tcell_subcluster.celltypes.umap.png"), width=13, height=6)

# subset it

subset_seu <- subset(seu_obj, subset = !tcell_clusters %in% c("6","7","11"))


DimPlot(subset_seu, reduction="umap.tcell_pca", group.by= c("tcell_clusters","cell_type"),
        label=T, label.box=T) 
ggsave(paste0(out_dir, "tcell_subcluster_filtered.celltypes.umap.png"), width=13, height=6)


SaveSeuratRds(subset_seu, file=paste0(out_dir, "tcell_subcluster.filtered.RDS"))
saveRDS(subset_seu)

# move annotation files to this directory
html_reports <- list.files(paste0(root_dir, "results/immune_cell_split/cassia_annotate_tcell/"),
                           "_report.html", full.names = T)
file.copy(from=html_reports, to=out_dir, overwrite = T)


