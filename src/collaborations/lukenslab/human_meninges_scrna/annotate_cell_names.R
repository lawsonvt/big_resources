library(Seurat)
library(SeuratData)
library(openxlsx)
library(ggplot2)

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/annotate_cell_names/")
dir.create(out_dir, showWarnings = F)

# read in seurat object
seu_obj <- LoadSeuratRds(paste0(root_dir, "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))

# read in annotations
cell_annot <- read.xlsx(paste0(root_dir, "cluster_names.final_8_24_2026.xlsx"))

# simplify cluster names
cell_annot$harmony_clusters <- as.character(as.numeric(gsub("cluster", "", cell_annot$harmony_clusters)))

# pull out metadata
metadata <- seu_obj@meta.data

cell_annot$harmony_clusters <- factor(as.character(cell_annot$harmony_clusters),
                                      levels = levels(metadata$harmony_clusters))

# merge in annotations
metadata <- merge(metadata,
                  cell_annot,
                  by="harmony_clusters")

# fix order
rownames(metadata) <- metadata$cell_id
metadata <- metadata[colnames(seu_obj),]

# add into seurat object
seu_obj@meta.data$claude_prediction <- metadata$claude_prediction
seu_obj@meta.data$gemini_prediction <- metadata$gemini_prediction
seu_obj@meta.data$final_cell_name <- metadata$final_names
seu_obj@meta.data$cell_category <- metadata$cell_category

# make some plots with names
DimPlot(seu_obj, reduction = "umap.harmony", group.by="harmony_clusters",
        label=T) + theme(legend.position = "none")

DimPlot(seu_obj, reduction = "umap.harmony", group.by="final_cell_name",
        label=T, label.size=3, label.box = T, repel=T)  + theme(legend.position = "none")
ggsave(paste0(out_dir, "final_cell_names.umap.png"), width=9, height=7)

DimPlot(seu_obj, reduction = "umap.harmony", group.by="cell_category",
        label=T, label.size=3, label.box = T, repel=T)  + theme(legend.position = "none")
ggsave(paste0(out_dir, "cell_category.umap.png"), width=9, height=7)

# cell count bar plots

ggplot(metadata,
       aes(x=condition,
           fill=Sample_name)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ final_names, ncol=5, scales="free_y") +
  labs(x=NULL, y="Cell Count", fill=NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "final_cell_names.bar_plot.png"), width=14, height=11)
  
# cell category bar plots

ggplot(metadata,
       aes(x=condition,
           fill=Sample_name)) +
  geom_bar(color="black") +
  theme_bw() +
  facet_wrap(~ cell_category, ncol=4, scales="free_y") +
  labs(x=NULL, y="Cell Count", fill=NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(out_dir, "cell_category.bar_plot.png"), width=12, height=9)

# save seurat object
SaveSeuratRds(seu_obj, paste0(out_dir, "seurat.named.RDS"))


