library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggsci)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/cell_clustering_plots/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

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

# count up results
cell_counts <- as.data.frame(table(metadata$condition,
                     metadata$cell_type))
colnames(cell_counts) <- c("condition", "cell_type", "count")

total_counts <- as.data.frame(table(metadata$condition))
colnames(total_counts) <- c("condition", "total_count")

cell_counts <- merge(cell_counts, total_counts,
                     by="condition")
cell_counts$fraction <- cell_counts$count / cell_counts$total_count

cell_counts <- cell_counts[order(cell_counts$fraction, decreasing=T),]
cell_counts$cell_type <- factor(as.character(cell_counts$cell_type),
                                levels=unique(as.character(cell_counts$cell_type)))

ggplot(cell_counts, 
       aes(x=condition,
           y=fraction,
           fill=cell_type)) +
  geom_bar(stat="identity", color="white") +
  theme_bw() +
  labs(x=NULL, y="Cell Count Proportion", fill=NULL)
ggsave(paste0(out_dir, "proportion_cell_types.stacked_bar_plot.png"), width=7, height=5)

cell_counts$cell_type <- factor(as.character(cell_counts$cell_type),
                                levels=rev(unique(as.character(cell_counts$cell_type))))

ggplot(cell_counts,
       aes(y=cell_type,
           x=fraction,
           fill=condition)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  labs(x="Cell Count Proportion", y=NULL, fill=NULL)
ggsave(paste0(out_dir, "proportion_cell_types.dodged_celltype_bar_plot.png"), width=7, height=5)

cell_counts$cell_type <- factor(as.character(cell_counts$cell_type),
                                levels=unique(as.character(cell_counts$cell_type)))

ggplot(cell_counts,
       aes(x=condition,
           y=count)) +
  geom_bar(stat="identity", color="black", fill="grey") +
  facet_wrap(~ cell_type, scale="free_y", ncol=6) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "cell_type_counts.facet_bar_plot.png"), width=12, height=8)
  
