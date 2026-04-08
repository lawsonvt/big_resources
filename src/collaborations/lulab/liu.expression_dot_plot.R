library(monocle3)
library(ggplot2)

out_dir <- "~/Documents/projects/lulab/yenchen_liu/plots/"

dir.create(out_dir, showWarnings = F)

# same file from analysis i did for casey
cds <- readRDS("~/Documents/projects/lulab/casey_bauchle/GSE249268_Final_aging_OPC_cds_mm106.rds")

metadata <- as.data.frame(colData(cds))
genedata <- as.data.frame(rowData(cds))


metadata_subset <- metadata[metadata$final_identity == "quiescent_OPC",]

# subset to quiescent OPC
cds_subset <- cds[, rownames(metadata_subset)]


# gene groups
p2x_genes <- c("P2rx1", "P2rx2", "P2rx3", "P2rx4", "P2rx5", "P2rx6", "P2rx7")
p2y_genes <- c("P2ry1", "P2ry2", "P2ry4", "P2ry6", "P2ry11", "P2ry12", "P2ry13", "P2ry14")

plot_genes_by_group(cds_subset,
                    p2x_genes,
                    group_cells_by = "new_age",
                    ordering_type = "none") +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "p2x_genes.gene_expression_dot_plot.png"),
       width=6, height=5)


plot_genes_by_group(cds_subset,
                    p2y_genes,
                    group_cells_by = "new_age",
                    ordering_type = "none") +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "p2y_genes.gene_expression_dot_plot.png"),
       width=6, height=5)

plot_genes_by_group(cds_subset,
                    c(p2x_genes,
                      p2y_genes),
                    group_cells_by = "new_age",
                    ordering_type = "none") +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "p2x_p2y_genes.gene_expression_dot_plot.png"),
       width=6, height=5)
ggsave(paste0(out_dir, "p2x_p2y_genes.gene_expression_dot_plot.svg"),
       width=6, height=5)
