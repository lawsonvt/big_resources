library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/integrate_subclustering_results/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))


# get subclustering directories
sub_dirs <- list.dirs(paste0(root_dir, "results/"), full.names=T, recursive = F)
sub_dirs <- grep("\\/subclustering", sub_dirs, value=T)

# doublet data
doublet_metas <- lapply(sub_dirs, function(d) {
  
  doublet_file <- list.files(d, ".doublet_output.RDS", full.names = T)
  
  data <- readRDS(doublet_file)
  data$cell_id <- rownames(data)
  
  return(data)
})
doublet_metas_df <- bind_rows(doublet_metas)

# subclustering metadata
subcluster_metas <- lapply(sub_dirs, function(d) {
  
  meta_file <- list.files(d, ".metadata.RDS", full.names = T)

  data <- readRDS(meta_file)
  data$cell_id <- rownames(data)
  
  return(data)
  
})
subcluster_metas_df <- bind_rows(subcluster_metas)

# get doublets and remove from seurat object
doublet_cell_ids <- doublet_metas_df[doublet_metas_df$scDblFinder.class == "doublet",]$cell_id

meta <- seu_obj@meta.data

total_cell_ids <- rownames(meta)

# remove doublet cell ids
total_cell_ids <- total_cell_ids[!total_cell_ids %in% doublet_cell_ids]

seu_obj <- subset(seu_obj, cells=total_cell_ids)

# new meta file
meta <- seu_obj@meta.data

# get data that was not part of subclustering
meta_no_subs <- meta[!meta$cell_type %in% subcluster_metas_df$cell_type,]

# now combine with the meta from subclustering (that was filtered)
total_cell_ids <- c(rownames(meta_no_subs),
                    subcluster_metas_df$cell_id)

# subset again, filtering out cells that seemed to be sample specific
seu_obj <- subset(seu_obj, cells=total_cell_ids)

# redo UMAP plots
DimPlot(seu_obj,
        reduction = "umap.harmony_filtered",
        group.by="cell_type", label=T) +
  labs(x="UMAP1", y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "cell_type.umap.total.png"), width=14, height=9)

# group the conditions
metadata <- seu_obj@meta.data

metadata$sample <- factor(metadata$orig.ident,
                          levels=c("X1161_1NEG",
                                   "X1176-WT1",
                                   "X1177-WT2",
                                   "X1162_2Het",
                                   "X1178-KO1",
                                   "X1179-KO2"))

metadata$condition <- "WT"

metadata[metadata$sample %in% c("X1162_2Het",
                                "X1178-KO1",
                                "X1179-KO2"),]$condition <- "KO"
metadata$condition <- factor(metadata$condition,
                             levels=c("WT","KO"))

# add to seurat object
seu_obj$sample <- metadata$sample
seu_obj$condition <- metadata$condition

# now make the plots per condition

DimPlot(seu_obj,
        reduction = "umap.harmony_filtered",
        group.by="cell_type", label=T, 
        label.size = 3,
        split.by="condition",
        ncol=2) +
  labs(x="UMAP1", y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "cell_type.umap.per_condition.png"), width=18, height=8)

# save seurat file
SaveSeuratRds(seu_obj, paste0(out_dir, "total_samples.post_subclustering.seurat.RDS"))




