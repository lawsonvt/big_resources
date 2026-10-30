library(Seurat)
library(ggplot2)

root_dir <- "~/Documents/projects/harrislab/x1247-GEX/"

# test comment

out_dir <- paste0(root_dir, "results/hto_demultiplex/")
dir.create(out_dir, showWarnings = F)

# read in the h5 file
h5_file <- paste0(root_dir, "cellranger_outputs.multi_cloud10x/filtered_feature_bc_matrix.h5")

raw_data <- Read10X_h5(h5_file)
# it contains BOTH!!!!!

# create the seurat object
seu <- CreateSeuratObject(counts=raw_data$`Gene Expression`, 
                          project="x1247-GEX", 
                          min.cells = 3, 
                          min.features = 200)

hto_data <- raw_data$`Antibody Capture`
common_cells <- intersect(colnames(seu), colnames(hto_data))

# subset to common cells

seu <- subset(seu, cells=common_cells)
hto_data <- hto_data[,common_cells]


seu[["HTO"]] <- CreateAssayObject(counts = hto_data)

seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")

seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)

metadata <- seu@meta.data

table(metadata$hash.ID, metadata$HTO_classification.global)

metadata$sample_id <- gsub("-TotalSeqB", "", metadata$hash.ID)

ggplot(metadata,
       aes(x=sample_id)) +
  geom_bar() +
  theme_bw()

# let's make this prettier
sample_counts <- as.data.frame(table(metadata$sample_id))
colnames(sample_counts) <- c("sample_id","cell_count")

sample_counts <- sample_counts[order(sample_counts$cell_count,
                                     decreasing=T),]
sample_counts$sample_id <- factor(as.character(sample_counts$sample_id),
                                  levels=as.character(sample_counts$sample_id))

ggplot(sample_counts,
       aes(x=sample_id,
           y=cell_count)) +
  geom_bar(stat="identity", color="black", fill="maroon") +
  geom_text(aes(label=cell_count), vjust=-0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "raw_sample_cell_counts.png"), width=7, height=5)

# save sample IDs to seurart

seu@meta.data$sample_id <- metadata$sample_id

# save Seurat
SaveSeuratRds(seu, file=paste0(out_dir, "seurat.hto_demultiplexed.RDS"))



