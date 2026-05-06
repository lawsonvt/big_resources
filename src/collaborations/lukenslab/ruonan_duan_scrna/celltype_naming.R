library(Seurat)
library(openxlsx)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/celltype_naming/")
dir.create(out_dir, showWarnings = F)


# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# cell annotations
cell_annot <- read.xlsx(paste0(root_dir, "cell type 20260329.xlsx"))
colnames(cell_annot) <- c("harmony_clusters","cell_type")

metadata <- int_seu@meta.data

metadata$cell_id <- rownames(metadata)

cell_annot$harmony_clusters <- factor(cell_annot$harmony_clusters,
                                     levels=levels(metadata$harmony_clusters))

# merge it in
metadata <- merge(metadata,
                  cell_annot,
                  by="harmony_clusters")

rownames(metadata) <- metadata$cell_id

metadata <- metadata[colnames(int_seu),]

int_seu$cell_type <- metadata$cell_type

# plot it out
# DimPlot(int_seu, reduction = "umap.harmony",
#         group.by = "harmony_clusters",
#         label=T)
# ggsave(paste0(out_dir, "harmony_clusters_v_cell_types.umap.png"), width=16, height=8)


DimPlot(int_seu, reduction = "umap.harmony",
                 group.by = "harmony_clusters",
                 label=T) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "harmony_clusters.umap.png"), width=8, height=6)


DimPlot(int_seu, reduction = "umap.harmony",
        group.by = "cell_type",
        label=T) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "celltype_clusters.umap.png"), width=12, height=8)

DimPlot(int_seu, reduction = "umap.harmony",
        group.by = "cell_type",
        split.by="orig.ident",
        label=F, ncol=2,
        raster=F) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "celltype_clusters.umap.by_sample.png"), width=14, height=12)


# drop a cluster

int_seu <- subset(int_seu, subset = harmony_clusters != "24")

# redo PCA/UMAP without dropped cluster
int_seu <- RunPCA(int_seu)

max_pc_dim <- 20

int_seu <- FindNeighbors(int_seu, dims=1:max_pc_dim)

# create new ummap
int_seu <- RunUMAP(int_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.harmony_filtered")

DimPlot(int_seu, reduction = "umap.harmony_filtered",
        group.by = "harmony_clusters",
        label=T) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "harmony_clusters.filtered_umap.png"), width=8, height=6)


DimPlot(int_seu, reduction = "umap.harmony_filtered",
        group.by = "cell_type",
        label=T) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "celltype_clusters.filtered_umap.png"), width=12, height=8)


DimPlot(int_seu, reduction = "umap.harmony_filtered",
        group.by = "cell_type",
        split.by="orig.ident",
        label=F, ncol=2,
        raster=F) +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "celltype_clusters.filtered_umap.by_sample.png"), width=14, height=12)

# save seurat object

SaveSeuratRds(int_seu, paste0(out_dir, "all_samples.celltype_named.seurat.RDS"))


# more plots
# group the conditions
metadata <- int_seu@meta.data

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
int_seu$sample <- metadata$sample
int_seu$condition <- metadata$condition


DimPlot(int_seu, reduction = "umap.harmony_filtered",
        group.by = "cell_type",
        label=T, split.by="condition") +
  labs(x="UMAP1",y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "celltype_clusters.filtered_umap.by_condition.png"), width=16, height=8)






