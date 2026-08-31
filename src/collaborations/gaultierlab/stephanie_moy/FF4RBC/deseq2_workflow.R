library(DESeq2)
library(openxlsx)
library(ggplot2)
library(gtools)
library(ggrepel)
library(sva)


root_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/"

out_dir <- paste0(root_dir, "FF4RBC/results/salmon_deg_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

exp_matrix <- read.delim(paste0(root_dir, "FF4RBC_results/FF4RBC-expression-matrix.tsv"))

metadata <- read.xlsx(paste0(root_dir, "FF4RBC/metadata.xlsx"))

# reformat counts
rownames(exp_matrix) <- exp_matrix$gene_id
any(duplicated(exp_matrix$gene_name))

gene_xref <- unique(exp_matrix[,c("gene_id","gene_name","gene_biotype")])

# if no gene name, just use ensembl id
gene_xref[gene_xref$gene_name == "",]$gene_name <-
  gene_xref[gene_xref$gene_name == "",]$gene_id


# counts
counts_mat <- exp_matrix[,mixedsort(grep("count", colnames(exp_matrix), value=T))]
cpm_mat <- exp_matrix[,mixedsort(grep("cpm", colnames(exp_matrix), value=T))]

# fix column names
colnames(counts_mat) <- gsub("_count","",colnames(counts_mat))

# convert counts to integers
counts_mat <- round(counts_mat)

# fix metadata
rownames(metadata) <- metadata$sample_id

metadata$condition <- factor(metadata$condition)
metadata$sex <- factor(metadata$sex)


# ensure same order
metadata <- metadata[colnames(counts_mat),]

# make DDS object
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = metadata,
                              design = ~condition)


# pre filter counts data (for plotting purposes)
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)


# QC plots
vsd <- vst(dds)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape=sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=sample_id), color="black") +
  theme_bw()
ggsave(paste0(out_dir, "pca_plot.png"), width=7, height=5)

# Get normalized counts for SVA
dds_norm <- estimateSizeFactors(dds)
norm_counts <- counts(dds_norm, normalized = TRUE)


mod <- model.matrix(~ 0 + condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))

# Calculate surrogate variables
svobj <- svaseq(norm_counts, mod, mod0)

# PCA after correction
vsd_corrected <- vst(dds, blind = FALSE)
assay(vsd_corrected) <- limma::removeBatchEffect(
  assay(vsd_corrected),
  covariates = svobj$sv,
  design = mod
)

pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))

ggplot(pcaData_after, aes(x = PC1, y = PC2, color = condition, shape=sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA - SVA applied") +
  geom_text_repel(aes(label=sample_id), color="black") +
  theme_bw()
ggsave(paste0(out_dir, "pca_plot.sva.png"), width=7, height=5)

# Add SVs to colData
for (i in 1:ncol(svobj$sv)) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

# Create new design formula including SVs
# Build the SV terms dynamically
sv_terms <- paste0("SV", 1:ncol(svobj$sv), collapse = " + ")
design_formula <- as.formula(paste("~", sv_terms, " + condition"))

# Update the design
design(dds) <- design_formula

dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","treatment","control"))

res <- as.data.frame(res)
# remove NAs
res <- res[!is.na(res$padj),]

# merge in gene names
res$gene_id <- rownames(res)

res <- merge(gene_xref,
             res, by="gene_id")

write.xlsx(res, paste0(out_dir, "deg_results.xlsx"), colWidths="auto")


# lets make a volcan0 plot!

subset <- res

subset$log_p <- -log10(subset$pvalue)

subset_sig <- subset[subset$padj < 0.05 &
                       abs(subset$log2FoldChange) > 0.5,]

logp_thresh <- min(subset_sig$log_p)

ggplot(subset,
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  geom_hline(yintercept = logp_thresh,
             color="red", linetype=2) +
  geom_vline(xintercept = 0.5,
             color="red", linetype=2) +
  geom_vline(xintercept = -0.5,
             color="red", linetype=2) +
  geom_point(data=subset_sig,
             color="red",
             alpha=0.4) +
  geom_text_repel(data=subset_sig,
                  aes(label=gene_name),
                  color="red", size=2.5,
                  max.overlaps = 50) +
  theme_bw() +
  labs(x="Log2 Fold Change", y="-log10(P-Value)")
ggsave(paste0(out_dir, "volcano_plot.png"), width=6, height=5)


