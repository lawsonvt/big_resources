library(DESeq2)
library(ggplot2)
library(readxl)
library(ggrepel)
library(sva)
library(limma)
library(tibble)
library(reshape2)
library(openxlsx)

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

# alternate PCA to find loadings
vsd_mat <- assay(vsd)

# Optional: Select top 500 most variable genes (DESeq2's default behavior)
rv <- rowVars(vsd_mat)
select <- order(rv, decreasing = TRUE)[seq_len(500)]
vsd_subset <- vsd_mat[select, ]

pca_vals <- prcomp(t(vsd_subset))

pc1_loadings <- sort(pca_vals$rotation[,1], decreasing=T)


# try out sva

# Get normalized counts for SVA
dds_norm <- estimateSizeFactors(dds)
norm_counts <- counts(dds_norm, normalized = TRUE)

# dat <- vst(dds_norm, blind = FALSE)
# dat_matrix <- assay(dat)

mod <- model.matrix(~ 0 + condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))

n_sv <- num.sv(norm_counts, mod)

# Calculate surrogate variables
svobj <- svaseq(norm_counts, mod, mod0, n.sv = n_sv)


vsd_corrected <- vst(dds, blind = FALSE)
assay(vsd_corrected) <- limma::removeBatchEffect(
  assay(vsd_corrected),
  covariates = svobj$sv,
  design = mod
)

pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))

ggplot(pcaData_after, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA, SVA applied") +
  geom_text_repel(aes(label=Name)) +
  theme_bw()
ggsave(paste0(out_dir, "pca_sva_plot.png"), width=6, height=5)

# try doing differential expression analysis

# Add SVs to colData
for (i in 1:n_sv) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

# Create new design formula including SVs
# Build the SV terms dynamically
sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")

design_formula <- as.formula(paste("~", sv_terms, "+ condition"))

# Update the design
design(dds) <- design_formula # no noticeable differences in results, so dont use

dds <- DESeq(dds)

# create the contrasts
contrasts <- list("FUS-Control"=c("condition","FUS","Control"),
                  "CAR-Control"=c("condition","CAR","Control"),
                  "CAR_FUS-Control"=c("condition","CAR_FUS","Control"))

results_list <- lapply(contrasts, function(contrast) {
  
  print(paste0(contrast[2], "-",
               contrast[3]))
  
  # calculate results
  res <- results(dds, contrast) 
  
  res <- as.data.frame(res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  
  # merge in gene names
  res <- rownames_to_column(res, var = "gene")
  
  # sort the data
  res <- res[order(res$pvalue),]
  
  res$contrast <- paste0(contrast[2], "-",
                         contrast[3])
  
  return(res)
  
})
names(results_list) <- names(contrasts)

write.xlsx(results_list, paste0(out_dir, "deg_results.xlsx"), colWidths="auto")

# make some plots
vsd_counts <- assay(vsd_corrected)

vsd_counts_long <- melt(vsd_counts)
colnames(vsd_counts_long) <- c("gene","sample_id","value")

vsd_counts_long$gene <- as.character(vsd_counts_long$gene)


vsd_counts_long <- merge(vsd_counts_long,
                         metadata, by="sample_id")

exp_genes <- unique(vsd_counts_long$gene)

grep("Mcp", exp_genes, value=T)

genelist <- c("Hspa4","Ifng", "Tnf","Il1b","Ccl2","Hgf")

# try just ong big plot

ggplot(vsd_counts_long[vsd_counts_long$gene %in% genelist,],
       aes(x=condition,
           y=value,
           color=condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(size=3) +
  theme_bw() +
  facet_wrap(~ gene, ncol=3, scales="free_y") +
  geom_text_repel(aes(label=Name), size=3) +
  guides(color="none") +
  labs(x=NULL, y="Normalized Expression")
ggsave(paste0(out_dir, "example_gene_expression_plots.png"), width=12, height=8)

# do the same for PC1 loadings values

ggplot(vsd_counts_long[vsd_counts_long$gene %in% names(head(pc1_loadings)),],
       aes(x=condition,
           y=value,
           color=condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(size=3) +
  theme_bw() +
  facet_wrap(~ gene, ncol=3, scales="free_y") +
  geom_text_repel(aes(label=Name), size=3) +
  guides(color="none") +
  labs(x=NULL, y="Normalized Expression")
ggsave(paste0(out_dir, "top_pc1_gene_expression_plots.png"), width=12, height=8)






