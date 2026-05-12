library(DESeq2)
library(ggplot2)
library(readxl)
library(ggrepel)

root_dir <- "~/Documents/projects/sheybanilab/stef_maslova/"

out_dir <- paste0(root_dir, "results/explore_counts/")
dir.create(out_dir, showWarnings = F, recursive = T)


# read in data
metadata <- as.data.frame(read_excel(paste0(root_dir, "260130_Bulk_Seq_List.xls")))
colnames(metadata)[1] <- "sample_id"

metadata$sample_id <- paste0("Mouse_", metadata$sample_id)
rownames(metadata) <- metadata$sample_id

metadata$condition <- factor(metadata$condition,
                             levels=unique(metadata$condition))


counts <- read.delim(paste0(root_dir, "D7_RNA_Seq/countsAllStef.txt"), comment.char = "#")

# clean up counts
rownames(counts) <- counts$Geneid

# split it up
gene_xref <- counts[,1:6]

counts_mat <- counts[,7:18]

# fix column names
colnames(counts_mat) <- sapply(colnames(counts_mat), function(x) {
  unlist(strsplit(x, "\\."))[6]
})

all(colnames(counts_mat) %in% rownames(metadata)) 
counts_mat <- counts_mat[,rownames(metadata)]

# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts_mat >= 10) >= 3
counts_mat <- counts_mat[keep, ]

# create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = metadata,
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# QC plots
vsd <- vst(dds)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=Name)) +
  theme_bw()
ggsave(paste0(out_dir, "pca_plot.png"), width=6, height=5)

# distance matrix
sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)


