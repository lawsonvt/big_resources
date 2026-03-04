library(DESeq2)
library(tximport)
library(ggplot2)
library(openxlsx)
library(snakecase)
library(data.table)
library(ggrepel)
library(cowplot)
library(sva)
library(limma)
library(pbayes)
library(tibble)

out_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/bulkRNASeq/results/deg_workflow/"

dir.create(out_dir, recursive = T, showWarnings = F)

# get metadata
metadata <- read.xlsx("~/Documents/projects/gaultierlab/stephanie_moy/bulkRNASeq/Sequencing info 25027.v2.xlsx")
metadata$sample_id <- paste0("Sample_", metadata$experiment.25027)
metadata$condition <- gsub("\\+", "", gsub(" ", "_", metadata$name))
metadata$condition <- factor(as.character(metadata$condition),
                             levels=c("WT_OPC",
                                      "KO_OPC",
                                      "WT_OL_T3",
                                      "KO_OL_T3"))
rownames(metadata) <- metadata$sample_id


# get the quant files
quant_files <- list.files("~/Documents/projects/gaultierlab/stephanie_moy/bulkRNASeq/quants/",
                          "quant.sf",
                          recursive=T, full.names = T)

names(quant_files) <- sapply(quant_files, function(x) {
  
  dirs <- unlist(strsplit(dirname(x), "\\/"))
  
  return(gsub("quant", "", last(dirs)))
  
})

# ensure that files and metadata have same order
all(names(quant_files) %in% rownames(metadata))
all(rownames(metadata) %in% names(quant_files) )

metadata <- metadata[names(quant_files),]

# get tx2gene
tx2gene <- read.csv("~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/Mus_musculus.GRCm39.transcript2gene.csv")

# read in quant files
txi <- tximport(quant_files, type="salmon", 
                tx2gene = tx2gene[,c("transcript_id","gene_id")])

# create DESeq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~0 + condition)

# pre filter counts data (for plotting purposes)
smallestGroupSize <- 4
keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
ddsTxi <- ddsTxi[keep,]

ddsTxi <- DESeq(ddsTxi)

# QC plots
vsd <- vst(ddsTxi)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()
ggsave(paste0(out_dir, "pca_plot.png"), width=6, height=5)


p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = condition)) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p2 <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = factor(batch))) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(fill="batch") +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p3 <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = Mbp)) +
  geom_point(size = 3, pch=21) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  scale_fill_gradient(low="white",high="purple") +
  theme_bw()

plot_grid(p1,p2,p3, ncol=2)
ggsave(paste0(out_dir, "pca_multiplot.png"), width=10, height=8, bg="white")

mod <- model.matrix(~ 0 + condition, colData(ddsTxi))

vsd_corrected <- vst(ddsTxi, blind = FALSE)
assay(vsd_corrected) <- limma::removeBatchEffect(
  assay(vsd_corrected),
  batch = factor(metadata$batch),
  design = mod
)

pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))

ggplot(pcaData_after, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p1 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = condition)) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p2 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = factor(batch))) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  labs(fill="batch") +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p3 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = Mbp)) +
  geom_point(size = 3, pch=21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  scale_fill_gradient(low="white",high="purple") +
  theme_bw()

plot_grid(p1,p2,p3, ncol=2)
ggsave(paste0(out_dir, "pca_multiplot.batch_removed.png"), width=10, height=8, bg="white")

# try SVA

# Get normalized counts for SVA
dds_norm <- estimateSizeFactors(ddsTxi)
norm_counts <- counts(dds_norm, normalized = TRUE)

dat <- vst(dds_norm, blind = FALSE)
dat_matrix <- assay(dat)

mod <- model.matrix(~ 0 + condition, colData(ddsTxi))
mod0 <- model.matrix(~1, colData(ddsTxi))

n_sv <- num.sv(dat_matrix, mod)

# Calculate surrogate variables
svobj <- sva(dat_matrix, mod, mod0, n.sv = n_sv)

# PCA after correction
vsd_corrected <- vst(ddsTxi, blind = FALSE)
assay(vsd_corrected) <- limma::removeBatchEffect(
  assay(vsd_corrected),
  covariates = svobj$sv,
  design = mod
)

pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))

p1 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = condition)) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA - SVA applied") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p2 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = factor(batch))) +
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  labs(fill="batch") +
  ggtitle("PCA - SVA applied") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  theme_bw()

p3 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, fill = Mbp)) +
  geom_point(size = 3, pch=21) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA - SVA applied") +
  geom_text_repel(aes(label=experiment.25027), color="black") +
  scale_fill_gradient(low="white",high="purple") +
  theme_bw()

plot_grid(p1,p2,p3, ncol=2)
ggsave(paste0(out_dir, "pca_multiplot.sva_applied.png"), width=10, height=8, bg="white")


# Add SVs to colData
for (i in 1:n_sv) {
  colData(ddsTxi)[, paste0("SV", i)] <- svobj$sv[, i]
}

# Create new design formula including SVs
# Build the SV terms dynamically
sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")
design_formula <- as.formula(paste("~", sv_terms, "+ condition"))

# Update the design
design(ddsTxi) <- design_formula

ddsTxi <- DESeq(ddsTxi)


# create the contrasts
contrasts <- list("KO_OPC-WT_OPC"=c("condition","KO_OPC","WT_OPC"),
                  "KO_OL_T3-WT_OL_T3"=c("condition","KO_OL_T3","WT_OL_T3"),
                  "KO_OPC-KO_OL_T3"=c("condition","KO_OPC","KO_OL_T3"),
                  "WT_OPC-WT_OL_T3"=c("condition","WT_OPC","WT_OL_T3"))

results_list <- lapply(contrasts, function(contrast) {
  
  # calculate results
  res <- results(ddsTxi, contrast) 
  
  res <- as.data.frame(res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  
  # posterior probability
  res$post_p <- pbayes(res$pvalue, n_cores=2, 
                   opt_method="SANN")$posterior_prob
  
  # merge in gene names
  res <- rownames_to_column(res, var = "gene_id")
  
  res <- merge(unique(tx2gene[,c("gene_id","gene_symbol")]),
               res,
               by="gene_id", all.y=T)
  
  # sort the data
  res <- res[order(res$pvalue),]
  
  res$contrast <- paste0(contrast[2], "-",
                         contrast[3])
  
  return(res)
  
})

# save the results
saveRDS(results_list, file=paste0(out_dir, "deg_results.sva.list.RDS"))

volcano_dir <- paste0(out_dir, "volcano_plots/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(results_list), function(contrast) {
  
  
  subset <- results_list[[contrast]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.05 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  p1 <- ggplot(subset,
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
    geom_text_repel(data=subset_top,
                    aes(label=gene_symbol),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=contrast)
  p1
  ggsave(paste0(volcano_dir, to_snake_case(contrast), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})

plot_grid(plotlist = volcano_plot_list, nrow = 2)

ggsave(paste0(out_dir, "all_contrasts.volcanoes.png"), width=7, height=6, bg="white")

