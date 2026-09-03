library(DESeq2)
library(openxlsx)
library(ggplot2)
library(gtools)
library(ggrepel)
library(sva)
library(reshape2)
library(ggVennDiagram)

# Preprocessing ----

root_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/"

out_dir <- paste0(root_dir, "FF4RBC/results/deseq2_workflow_sva_tests.permutations/")
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


permutes <- 100

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
  
  # Total SVA (no sex) ------
  
  # make DDS object
  dds_total <- DESeqDataSetFromMatrix(countData = counts_mat,
                                      colData = p_meta,
                                      design = ~condition)
  
  # pre filter counts data
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds_total) >= 10) >= smallestGroupSize
  dds_total <- dds_total[keep,]
  
  dds_total <- DESeq(dds_total)
  
  # Get normalized counts for SVA
  dds_norm <- estimateSizeFactors(dds_total)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  
  
  mod <- model.matrix(~ condition, colData(dds_total))
  mod0 <- model.matrix(~1, colData(dds_total))
  
  # Calculate surrogate variables
  svobj <- svaseq(norm_counts, mod, mod0)
  
  # Add SVs to colData
  for (i in 1:ncol(svobj$sv)) {
    colData(dds_total)[, paste0("SV", i)] <- svobj$sv[, i]
  }
  
  # Build the SV terms dynamically
  sv_terms <- paste0("SV", 1:ncol(svobj$sv), collapse = " + ")
  design_formula <- as.formula(paste("~ ", sv_terms, " + condition"))
  
  # Update the design
  design(dds_total) <- design_formula
  
  dds_total <- DESeq(dds_total)
  
  res_total <- results(dds_total, contrast=c("condition","treatment","control"))
  
  res_total <- as.data.frame(res_total)
  # remove NAs
  res_total <- res_total[!is.na(res_total$padj),]
  
  # merge in gene names
  res_total$gene_id <- rownames(res_total)
  
  res_total <- merge(gene_xref,
                     res_total, by="gene_id")
  
  res_total <- res_total[order(res_total$pvalue),]
  
  # Free SVA (sex included, but no SV restriction) ------
  
  # make DDS object
  dds_free <- DESeqDataSetFromMatrix(countData = counts_mat,
                                     colData = p_meta,
                                     design = ~sex+condition)
  
  # pre filter counts data
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds_free) >= 10) >= smallestGroupSize
  dds_free <- dds_free[keep,]
  
  dds_free <- DESeq(dds_free)
  
  # Get normalized counts for SVA
  dds_norm <- estimateSizeFactors(dds_free)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  
  
  mod <- model.matrix(~ sex+condition, colData(dds_free))
  mod0 <- model.matrix(~sex, colData(dds_free))
  
  # Calculate surrogate variables
  svobj <- svaseq(norm_counts, mod, mod0)
  
  # Add SVs to colData
  for (i in 1:ncol(svobj$sv)) {
    colData(dds_free)[, paste0("SV", i)] <- svobj$sv[, i]
  }
  
  # Build the SV terms dynamically
  sv_terms <- paste0("SV", 1:ncol(svobj$sv), collapse = " + ")
  design_formula <- as.formula(paste("~ sex + ", sv_terms, " + condition"))
  
  # Update the design
  design(dds_free) <- design_formula
  
  dds_free <- DESeq(dds_free)
  
  res_free <- results(dds_free, contrast=c("condition","treatment","control"))
  
  res_free <- as.data.frame(res_free)
  # remove NAs
  res_free <- res_free[!is.na(res_free$padj),]
  
  # merge in gene names
  res_free$gene_id <- rownames(res_free)
  
  res_free <- merge(gene_xref,
                    res_free, by="gene_id")
  
  res_free <- res_free[order(res_free$pvalue),]
  
  # 1 SV SVA (sex included, SV restricted to 1) ------
  
  # make DDS object
  dds_1sv <- DESeqDataSetFromMatrix(countData = counts_mat,
                                    colData = p_meta,
                                    design = ~sex+condition)
  
  # pre filter counts data
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds_1sv) >= 10) >= smallestGroupSize
  dds_1sv <- dds_1sv[keep,]
  
  dds_1sv <- DESeq(dds_1sv)
  
  # Get normalized counts for SVA
  dds_norm <- estimateSizeFactors(dds_1sv)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  
  
  mod <- model.matrix(~ sex+condition, colData(dds_1sv))
  mod0 <- model.matrix(~sex, colData(dds_1sv))
  
  # Calculate surrogate variables
  svobj <- svaseq(norm_counts, mod, mod0, n.sv = 1)
  
  # Add SVs to colData
  for (i in 1:ncol(svobj$sv)) {
    colData(dds_1sv)[, paste0("SV", i)] <- svobj$sv[, i]
  }
  
  # Build the SV terms dynamically
  sv_terms <- paste0("SV", 1:ncol(svobj$sv), collapse = " + ")
  design_formula <- as.formula(paste("~ sex + ", sv_terms, " + condition"))
  
  # Update the design
  design(dds_1sv) <- design_formula
  
  dds_1sv <- DESeq(dds_1sv)
  
  res_1sv <- results(dds_1sv, contrast=c("condition","treatment","control"))
  
  res_1sv <- as.data.frame(res_1sv)
  # remove NAs
  res_1sv <- res_1sv[!is.na(res_1sv$padj),]
  
  # merge in gene names
  res_1sv$gene_id <- rownames(res_1sv)
  
  res_1sv <- merge(gene_xref,
                   res_1sv, by="gene_id")
  
  res_1sv <- res_1sv[order(res_1sv$pvalue),]
  
  # Overlap analysis --------
  
  deg_direct_list <- list(
    "Total SVA" = paste0(res_total[res_total$padj < 0.05,]$gene_id,"|",
                         sign(res_total[res_total$padj < 0.05,]$log2FoldChange)),
    "Free SVA" = paste0(res_free[res_free$padj < 0.05,]$gene_id,"|",
                        sign(res_free[res_free$padj < 0.05,]$log2FoldChange)),
    "1SV SVA" = paste0(res_1sv[res_1sv$padj < 0.05,]$gene_id,"|",
                       sign(res_1sv[res_1sv$padj < 0.05,]$log2FoldChange))
  )
  
  data.frame(
    permutation = i,
    total_deg = length(deg_direct_list$`Total SVA`),
    free_deg = length(deg_direct_list$`Free SVA`),
    "1sv_deg" = length(deg_direct_list$`1SV SVA`),
    overlap = sum(deg_direct_list$`Total SVA` %in% deg_direct_list$`Free SVA` &
                    deg_direct_list$`Total SVA` %in% deg_direct_list$`1SV SVA`)
  )
})
permute_stats <- bind_rows(permute_stats)

summary(permute_stats)

sum(24 > permute_stats$overlap)


