library(DESeq2)
library(ggplot2)
library(readxl)
library(ggrepel)
library(sva)
library(limma)
library(tibble)
library(reshape2)
library(openxlsx)
library(snakecase)
library(cowplot)

root_dir <- "~/Documents/projects/sheybanilab/stef_maslova/"

out_dir <- paste0(root_dir, "results/hisat2_deg_analysis.drop_20rr/")
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

# drop 20RR
metadata <- metadata[metadata$Name != "20RR FUS",]

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

# try out sva

# Get normalized counts for SVA
dds_norm <- estimateSizeFactors(dds)
norm_counts <- counts(dds_norm, normalized = TRUE)

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

# drop results with low counts
results_list <- lapply(results_list, function(data) {
  
  return(data[data$baseMean > 30,])
  
})

write.xlsx(results_list, paste0(out_dir, "deg_results.xlsx"), colWidths="auto")

# save results to RDS
saveRDS(results_list, file=paste0(out_dir, "deg_results.RDS"))


# Make Volcano plots -----------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

cap <- 4

volcano_plot_list <-  lapply(names(results_list), function(contrast) {
  
  print(contrast)
  
  subset <- results_list[[contrast]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  # cap foldchange
  if (any(subset$log2FoldChange > cap)) {
    subset[subset$log2FoldChange > cap,]$log2FoldChange <- cap 
  }
  
  if (any(subset$log2FoldChange < -cap)) {
    subset[subset$log2FoldChange < -cap,]$log2FoldChange <- -cap 
  }
  
  subset_sig <- subset[subset$padj < 0.1,]
  
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
    # geom_vline(xintercept = 0.5,
    #            color="red", linetype=2) +
    # geom_vline(xintercept = -0.5,
    #            color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=gene),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=contrast)
  p1
  ggsave(paste0(volcano_dir, to_snake_case(contrast), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, nrow = 1)
ggsave(paste0(out_dir, "contrast.volcanoes.png"), width=10, height=5, bg="white")


# make some plots
vsd_counts <- assay(vsd_corrected)

vsd_counts_long <- melt(vsd_counts)
colnames(vsd_counts_long) <- c("gene","sample_id","value")

vsd_counts_long$gene <- as.character(vsd_counts_long$gene)


vsd_counts_long <- merge(vsd_counts_long,
                         metadata, by="sample_id")

top_gene_list <- results_list$`CAR_FUS-Control`$gene[1:9]

ggplot(vsd_counts_long[vsd_counts_long$gene %in% top_gene_list,],
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
ggsave(paste0(out_dir, "top_fdr_plots.car_fus.png"), width=8, height=7)

# same plots for other comparisons
top_gene_list <- results_list$`FUS-Control`$gene[1:9]

ggplot(vsd_counts_long[vsd_counts_long$gene %in% top_gene_list,],
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
ggsave(paste0(out_dir, "top_fdr_plots.fus.png"), width=8, height=7)


top_gene_list <- results_list$`CAR-Control`$gene[1:3]

ggplot(vsd_counts_long[vsd_counts_long$gene %in% top_gene_list,],
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
ggsave(paste0(out_dir, "top_fdr_plots.car.png"), width=8, height=3)




