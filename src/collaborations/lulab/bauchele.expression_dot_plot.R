library(monocle3)
library(ggplot2)

out_dir <- "~/Documents/projects/lulab/casey_bauchle/plots/"

dir.create(out_dir, showWarnings = F)

cds <- readRDS("~/Documents/projects/lulab/casey_bauchle/GSE249268_Final_aging_OPC_cds_mm106.rds")

metadata <- as.data.frame(colData(cds))
genedata <- as.data.frame(rowData(cds))

genes <- c("Adgrb1", "Elmo1", "Dock1", "Rac1")

metadata_subset <- metadata[metadata$final_identity == "quiescent_OPC",]

# subset to quiescent OPC
cds_subset <- cds[, rownames(metadata_subset)]

plot_genes_by_group(cds_subset,
                    genes,
                    group_cells_by = "new_age",
                    ordering_type = "none") +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "bai1_pathway.gene_expression_dot_plot.png"),
       width=6, height=5)
ggsave(paste0(out_dir, "bai1_pathway.gene_expression_dot_plot.svg"),
       width=6, height=5)
