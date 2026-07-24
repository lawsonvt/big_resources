library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)
library(cowplot)
library(gtools)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/endothelial_sub_diff_exp_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)


# read in DEGs
endo_degs <- readRDS(paste0(root_dir, "results/endothelial_subcluster.diff_exp_analysis/eae_minus_adult.clusters.de_results.RDS"))


# group things for heatmaps
endo_degs_df <- bind_rows(endo_degs)

table(endo_degs_df[endo_degs_df$padj < 0.05,]$cluster)


# get groups

# exam the genes of interest
group1_genes <- c("Il1r1",
                  "Tnfrsf1a", 
                  "Ifngr1",
                  "Ifngr2")


group1_degs <- endo_degs_df[endo_degs_df$gene %in% group1_genes,]


fc_matrix <- acast(group1_degs, cluster ~ gene, value.var="log2FoldChange")
fdr_matrix <- acast(group1_degs, cluster ~ gene, value.var="padj")

fdr_print <- fdr_matrix
fdr_print[fdr_matrix < 0.05] <- "*"
fdr_print[fdr_matrix > 0.05] <- ""
fdr_print[is.na(fdr_matrix)] <- ""

h1 <- Heatmap(fc_matrix,
              col=colorRamp2(c(-2,0,2), c("blue","white","red")),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(fdr_print[i,j], x, y)
              }, name="log2FC")
pdf(paste0(out_dir, "group1_heatmap.pdf"), width=6, height=5)
h1
dev.off()


group2_genes <- c("B2m",
                  "Tap1",
                  "Cxcl10",
                  "H2-Ab1")


group2_degs <- endo_degs_df[endo_degs_df$gene %in% group2_genes,]


fc_matrix <- acast(group2_degs, cluster ~ gene, value.var="log2FoldChange")
fdr_matrix <- acast(group2_degs, cluster ~ gene, value.var="padj")

fdr_print <- fdr_matrix
fdr_print[fdr_matrix < 0.05] <- "*"
fdr_print[fdr_matrix > 0.05] <- ""
fdr_print[is.na(fdr_matrix)] <- ""

h1 <- Heatmap(fc_matrix,
              col=colorRamp2(c(-2,0,2), c("blue","white","red")),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(fdr_print[i,j], x, y)
              }, name="log2FC")
pdf(paste0(out_dir, "group2_heatmap.pdf"), width=6, height=5)
h1
dev.off()


# group 3

group3_genes <- c("Stat1","Ciita", "Icam1", "Vcam1")

group3_degs <- endo_degs_df[endo_degs_df$gene %in% group3_genes,]


fc_matrix <- acast(group3_degs, cluster ~ gene, value.var="log2FoldChange")
fdr_matrix <- acast(group3_degs, cluster ~ gene, value.var="padj")

fdr_print <- fdr_matrix
fdr_print[fdr_matrix < 0.05] <- "*"
fdr_print[fdr_matrix > 0.05] <- ""
fdr_print[is.na(fdr_matrix)] <- ""

h1 <- Heatmap(fc_matrix,
              col=colorRamp2(c(-2,0,2), c("blue","white","red")),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(fdr_print[i,j], x, y)
              }, name="log2FC")
pdf(paste0(out_dir, "group3_heatmap.pdf"), width=6, height=5)
h1
dev.off()






