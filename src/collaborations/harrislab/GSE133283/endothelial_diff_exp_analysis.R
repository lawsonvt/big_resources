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

out_dir <- paste0(root_dir, "results/endothelial_diff_exp_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in cassia results
cassia_results <- readRDS(paste0(root_dir, "results/cassia_results_analysis/cassia_results_comparison.RDS"))

# read in DEGs
deg_results <- readRDS(paste0(root_dir, "results/cluster.diff_exp_analysis/eae_minus_adult.clusters.de_results.RDS"))

# subset the endothelial cells
endo_cells <- cassia_results[grepl("[Ee]ndothelial", cassia_results$gemini_prediction),]

endo_degs <- deg_results[endo_cells$harmony_clusters]

# volcano plots

volcano_dir <- paste0(out_dir, "volcano_plots.clusters/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(endo_degs), function(cluster) {
  
  print(cluster)
  
  subset <- endo_degs[[cluster]]
  
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
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=gene),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=cluster,
         subtitle="EAE - Adult")
  p1
  ggsave(paste0(volcano_dir, to_snake_case(cluster), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})

plot_grid(plotlist = volcano_plot_list, nrow = 2)
ggsave(paste0(out_dir, "endothelial_cluster_deg.volcanoes.png"), width=10, height=7, bg="white")

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


# output results
write.xlsx(endo_degs, paste0(out_dir, "endothelial_cluster_degs.xlsx"),
           colWidths="auto")

# output the endothelial cells
wb <- createWorkbook()
addWorksheet(wb, "comparison")
writeData(wb, "comparison", endo_cells)

n_rows <- nrow(endo_cells) + 1  # +1 for header row

setColWidths(wb, "comparison", cols=1, width=20)
setColWidths(wb, "comparison", cols=2:ncol(endo_cells), width=40)

setRowHeights(wb, "comparison", rows = 2:n_rows, heights = 80)

wrap_style <- createStyle(wrapText = TRUE)
addStyle(wb, "comparison",
         style = wrap_style,
         rows = 2:n_rows,
         cols = 2:ncol(endo_cells),       # whichever columns should wrap
         gridExpand = TRUE)    # applies style to every row/col combination

saveWorkbook(wb, paste0(out_dir, "endothelial_clusters.xlsx"), overwrite = T)




