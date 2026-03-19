library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pbayes)
library(openxlsx)

out_dir <- "~/Documents/projects/gaultierlab/stephanie_moy/bulkRNASeq/results/joint_p_analysis/"

dir.create(out_dir, showWarnings = F)

deg_list <- readRDS("~/Documents/projects/gaultierlab/stephanie_moy/bulkRNASeq/results/deg_workflow/deg_results.sva.list.RDS")

# add in posterior probs
opc_degs <- deg_list[["KO_OPC-WT_OPC"]]

opc_degs$post_p <- pbayes(opc_degs$pvalue, n_cores=2, level_pvals = T)$posterior_prob

olt3_degs <- deg_list[["KO_OL_T3-WT_OL_T3"]]

olt3_degs$post_p <- pbayes(olt3_degs$pvalue, n_cores=2, level_pvals = T)$posterior_prob
# combine ko v wt

joint_deg <- merge(opc_degs,
                   olt3_degs,
                   by=c("gene_id","gene_symbol"),
                   all=T,
                   suffixes=c(".opc",".olt3"))

joint_deg$joint_p <- joint_deg$post_p.opc *
  joint_deg$post_p.olt3

joint_deg$same_dir <- sign(joint_deg$log2FoldChange.opc) ==
  sign(joint_deg$log2FoldChange.olt3)

# differences
joint_deg$opc_not_olt3 <- joint_deg$post_p.opc * (1 - joint_deg$post_p.olt3)
joint_deg$olt3_not_opc <- joint_deg$post_p.olt3 * (1 - joint_deg$post_p.opc)

# make some heatmaps

fc_matrix <- as.matrix(joint_deg[,c("log2FoldChange.opc",
                                    "log2FoldChange.olt3")])
rownames(fc_matrix) <- joint_deg$gene_symbol
colnames(fc_matrix) <- c("OPC", "OL +T3")

fdr_matrix <- as.matrix(joint_deg[,c("padj.opc",
                                    "padj.olt3")])
rownames(fdr_matrix) <- joint_deg$gene_symbol
colnames(fdr_matrix) <- c("OPC", "OL +T3")

fdr_matrix.print <- fdr_matrix
fdr_matrix.print[fdr_matrix < 0.05] <- "*"
fdr_matrix.print[fdr_matrix > 0.05] <- ""

# make some heatmaps
top_same <- joint_deg[joint_deg$joint_p > 0.7,]$gene_symbol

subset_fc <- fc_matrix[top_same,]
subset_fdr <- fdr_matrix.print[top_same,]

ht <- Heatmap(subset_fc,
        col=colorRamp2(c(-2,0,2), c("blue","white","red")),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(subset_fdr[i,j], x, y)
        }, name="log2FC",
        cluster_columns = F)

pdf(paste0(out_dir, "joint_degs.heatmap.pdf"), width=5, height=4)
draw(ht, padding = unit(c(2, 2, 10, 2), "mm"))  # extra bottom padding

# Add footnote text at the bottom
grid.text(
  "'*' Adjusted p-value < 0.05",
  x = 0.8,
  y = unit(4, "mm"),          # distance from bottom
  just = "bottom",
  gp = gpar(fontsize = 10, fontface = "italic")  # optional styling
)
dev.off()

# in OPC not in OL-T3

top_opc <- joint_deg[joint_deg$opc_not_olt3 > 0.9 &
                       abs(joint_deg$log2FoldChange.opc) > 0.5,]$gene_symbol

subset_fc <- fc_matrix[top_opc,]
subset_fdr <- fdr_matrix.print[top_opc,]

ht <- Heatmap(subset_fc,
              col=colorRamp2(c(-2,0,2), c("blue","white","red")),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(subset_fdr[i,j], x, y)
              }, name="log2FC",
              cluster_columns = F)

pdf(paste0(out_dir, "opc_unique_degs.heatmap.pdf"), width=5, height=4)
draw(ht, padding = unit(c(2, 2, 10, 2), "mm"))  # extra bottom padding

# Add footnote text at the bottom
grid.text(
  "'*' Adjusted p-value < 0.05",
  x = 0.8,
  y = unit(4, "mm"),          # distance from bottom
  just = "bottom",
  gp = gpar(fontsize = 10, fontface = "italic")  # optional styling
)
dev.off()

# OL +T3

top_olt3 <- joint_deg[joint_deg$olt3_not_opc > 0.9 &
                       abs(joint_deg$log2FoldChange.olt3) > 0.5,]$gene_symbol

subset_fc <- fc_matrix[top_olt3,]
subset_fdr <- fdr_matrix.print[top_olt3,]

ht <- Heatmap(subset_fc,
              col=colorRamp2(c(-2,0,2), c("blue","white","red")),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(subset_fdr[i,j], x, y)
              }, name="log2FC",
              cluster_columns = F)

pdf(paste0(out_dir, "olt3_unique_degs.heatmap.pdf"), width=5, height=4)
draw(ht, padding = unit(c(2, 2, 10, 2), "mm"))  # extra bottom padding

# Add footnote text at the bottom
grid.text(
  "'*' Adjusted p-value < 0.05",
  x = 0.8,
  y = unit(4, "mm"),          # distance from bottom
  just = "bottom",
  gp = gpar(fontsize = 10, fontface = "italic")  # optional styling
)
dev.off()

# output as one excel
outlist <- list(joint_deg=joint_deg[order(joint_deg$joint_p,
                                          decreasing=T),],
                opc_unique=joint_deg[order(joint_deg$opc_not_olt3,
                                           decreasing = T),],
                olt3_unique=joint_deg[order(joint_deg$olt3_not_opc,
                                            decreasing = T),])

write.xlsx(outlist, paste0(out_dir, "joint_p_deg_results.xlsx"),
           colWidths="auto")

saveRDS(outlist, paste0(out_dir, "joint_p_deg_results.RDS"))




