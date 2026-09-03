library(DESeq2)
library(openxlsx)
library(ggplot2)
library(gtools)
library(ggrepel)
library(sva)
library(reshape2)
library(dplyr)

root_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/"

out_dir <- paste0(root_dir, "FF4RBC/results/deseq2_workflow.permutation_test/")
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

permutes <- 200

set.seed(42)

# function for determining if full rank
is_full_rank <- function(sex, condition) {
  design <- model.matrix(~ sex + condition)
  qr(design)$rank == ncol(design)
}

permute_stats <- lapply(1:permutes, function(i) {
  
  print(i)
  
  p_meta <- metadata
  p_meta$condition <- sample(p_meta$condition)
  
  if (!is_full_rank(p_meta$sex,
               p_meta$condition)) {
    p_meta$condition <- sample(p_meta$condition)
  }
  
  # make DDS object
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData = p_meta,
                                design = ~sex+condition)
  
  # pre filter counts data (for plotting purposes)
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  
  # Get normalized counts for SVA
  dds_norm <- estimateSizeFactors(dds)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  
  
  mod <- model.matrix(~ sex + condition, colData(dds))
  mod0 <- model.matrix(~sex, colData(dds))
  
  # Calculate surrogate variables
  svobj <- svaseq(norm_counts, mod, mod0, n.sv=1)
  
  # Add SVs to colData
  for (i in 1:ncol(svobj$sv)) {
    colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
  }
  
  # Create new design formula including SVs
  # Build the SV terms dynamically
  sv_terms <- paste0("SV", 1:ncol(svobj$sv), collapse = " + ")
  design_formula <- as.formula(paste("~ sex +", sv_terms, " + condition"))
  
  # Update the design
  design(dds) <- design_formula
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast=c("condition","treatment","control"))
  
  res <- as.data.frame(res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  
  # output values
  data.frame(permutation=i,
             nsv=svobj$n.sv,
             deg_count=nrow(res[res$padj < 0.05,]))
  
  
})
permute_stats <- bind_rows(permute_stats)

unp_degs <- 78

median(permute_stats$deg_count) / unp_degs
mean(permute_stats$deg_count)

sum(permute_stats$deg_count < unp_degs)


