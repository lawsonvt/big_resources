library(Seurat)
library(SeuratObject)
library(ggplot2)
library(openxlsx)

root_dir <- "/Users/mjl3p/Documents/projects/gaultierlab/shreyas_gujar/snRNASeq/"
#root_dir <- "/home/mjl3p/projects/gaultierlab/shreyas_gujar/snRNASeq/"

out_dir <- paste0(root_dir, "results/cell_type_naming/")
dir.create(out_dir, showWarnings = F)

seu_obj <- readRDS(paste0(root_dir, "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

cell_annotations <- read.xlsx(paste0(root_dir, "results/cassia_annotate/celltype_labels.xlsx"))

metadata <- seu_obj@meta.data

metadata$cell_id <- rownames(metadata)

# merge in cell names
metadata <- merge(metadata,
                  cell_annotations,
                  by.x="harmony_clusters",
                  by.y="Cluster.ID")

# fix order
rownames(metadata) <- metadata$cell_id

metadata <- metadata[colnames(seu_obj),]

# add in a condition
metadata$condition <- gsub("[0-9]$", "", metadata$orig.ident)
metadata$sample_id <- metadata$orig.ident

seu_obj@meta.data <- metadata

# DimPlot(seu_obj, reduction="umap.harmony", group.by=c("seurat_clusters", "cell_name"),
#         label=T)
# ggsave(paste0(out_dir, "clusters.umap.png"), width=18, height=8)
# 
# DimPlot(seu_obj, reduction="umap", group.by=c("seurat_clusters", "cell_name"))
# ggsave(paste0(out_dir, "clusters.umap.no_labels.png"), width=16, height=6)



DimPlot(seu_obj, reduction="umap.harmony", group.by="cell_name", label=T,
        label.box = T, repel = T)
ggsave(paste0(out_dir, "cell_names.umap.png"), width=12, height=9)

DimPlot(seu_obj, reduction="umap.harmony", group.by="cell_name")
ggsave(paste0(out_dir, "cell_names.umap.no_labels.png"), width=12, height=9)

DimPlot(seu_obj, reduction="umap.harmony", group.by="harmony_clusters", label=T,
        label.box = T, repel = T)
ggsave(paste0(out_dir, "clusters.umap.png"), width=10, height=9)

DimPlot(seu_obj, reduction="umap.harmony", group.by="harmony_clusters")
ggsave(paste0(out_dir, "clusters.umap.no_labels.png"), width=10, height=9)

DimPlot(seu_obj, reduction="umap.harmony", group.by="cell_name",
        split.by="sample_id", ncol=3, raster=F)
ggsave(paste0(out_dir, "cell_names.umap.samples.png"), width=14, height=9)

DimPlot(seu_obj, reduction="umap.harmony", group.by="cell_name",
        split.by="condition", ncol=2, raster=F)
ggsave(paste0(out_dir, "cell_names.umap.conditions.png"), width=13, height=6)

SaveSeuratRds(seu_obj, paste0(out_dir, "seurat.cell_named.RDS"))

# save metadata (since seurat is so massivve)
saveRDS(metadata, paste0(out_dir, "cell_named.metadata.RDS"))

# make some bar plots
sample_counts <- as.data.frame(table(metadata$sample_id), stringsAsFactors = F)
colnames(sample_counts) <- c("sample_id","sample_count")

condition_counts <- as.data.frame(table(metadata$condition), stringsAsFactors = F)
colnames(condition_counts) <- c("condition","condition_count")

sample_cell_counts <- as.data.frame(table(metadata$sample_id,
                                          metadata$cell_name),
                                    stringsAsFactors = F)
colnames(sample_cell_counts) <- c("sample_id","cell_name","sample_cell_count")

condition_cell_counts <- as.data.frame(table(metadata$condition,
                                             metadata$cell_name),
                                       stringsAsFactors = F)
colnames(condition_cell_counts) <- c("condition","cell_name","condition_cell_count")

# merge in total counts
sample_cell_counts <- merge(sample_cell_counts,
                            sample_counts,
                            by="sample_id")
sample_cell_counts$sample_percent <- (sample_cell_counts$sample_cell_count / 
  sample_cell_counts$sample_count) * 100


condition_cell_counts <- merge(condition_cell_counts,
                            condition_counts,
                            by="condition")
condition_cell_counts$condition_percent <- (condition_cell_counts$condition_cell_count / 
                                        condition_cell_counts$condition_count) * 100


# PLOTS!
ggplot(sample_counts, 
       aes(x=sample_id,
           y=sample_count)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  geom_text(aes(label=sample_count), vjust=-0.3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "total_sample_counts.png"), width=5, height=4)


sample_counts$condition <- gsub("[0-9]$", "", sample_counts$sample_id)

ggplot(sample_counts, 
       aes(x=condition,
           y=sample_count,
           fill=sample_id)) +
  geom_bar(stat="identity", color="black") +
  #geom_text(aes(label=sample_count), vjust=-0.3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "total_condition_counts.sample_split.png"), width=4, height=4)

ggplot(condition_counts, 
       aes(x=condition,
           y=condition_count)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  geom_text(aes(label=condition_count), vjust=-0.3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "total_condition_counts.png"), width=4, height=4)


ggplot(sample_cell_counts, 
       aes(x=sample_id,
           y=sample_cell_count)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  facet_wrap(~ cell_name, scales="free_y", ncol=5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "per_cell_sample_counts.png"), width=14, height=10)


ggplot(condition_cell_counts, 
       aes(x=condition,
           y=condition_cell_count)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  facet_wrap(~ cell_name, scales="free_y", ncol=5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "per_cell_condition_counts.png"), width=14, height=10)

ggplot(sample_cell_counts, 
       aes(x=sample_id,
           y=sample_percent)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  facet_wrap(~ cell_name, scales="free_y", ncol=4) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Percent of Sample Cells")
ggsave(paste0(out_dir, "per_cell_sample_percent.png"), width=14, height=10)


ggplot(condition_cell_counts, 
       aes(x=condition,
           y=condition_percent)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  facet_wrap(~ cell_name, scales="free_y", ncol=4) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Percent of Condition Cells")
ggsave(paste0(out_dir, "per_cell_condition_percent.png"), width=14, height=10)

