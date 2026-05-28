library(monocle3)
library(ggplot2)
library(reshape2)

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

# another genelist

genes2 <- c("Lrp1",
            "Adgrb1",
            "Axl",
            "Mertk",
            "Megf10",
            "C3",
            "Trem2",
            "Adgrg1",
            "Cx3cr1")

plot_genes_by_group(cds_subset,
                    genes2,
                    group_cells_by = "new_age",
                    ordering_type = "none") +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "gene_list.gene_expression_dot_plot.png"),
       width=6, height=5)
ggsave(paste0(out_dir, "gene_list.gene_expression_dot_plot.svg"),
       width=6, height=5)

metadata_subset2 <- metadata_subset[metadata_subset$new_age == "1_P30",]
genedata_subset2 <- genedata[genedata$gene_short_name %in% genes2,]


cds_subset2 <- cds_subset[rownames(genedata_subset2), rownames(metadata_subset2)]



plot_genes_violin(cds_subset2,
                  group_cells_by = "new_age",
                  nrow=1, ncol=9) + 
  geom_jitter(alpha=0.2, size=0.1, width=0.2)

counts <- as.matrix(normalized_counts(cds_subset2, norm_method = "size_only"))

#counts <- SingleCellExperiment::counts(cds_subset2)
#counts <- Matrix::t(Matrix::t(counts)/size_factors(cds_subset2))

counts_long <- melt(as.matrix(counts))
colnames(counts_long) <- c("gene_id","barcode","value")

counts_long <- merge(counts_long,
                     genedata,
                     by.x="gene_id",
                     by.y="gene_id.x")

ggplot(counts_long,
       aes(x=gene_short_name,
           y=value)) +
  geom_violin(fill="lightblue") +
  theme_bw() + scale_y_log10() +
  geom_jitter(width=0.1, size=0.1, alpha=0.2) +
  labs(x=NULL, y="Expression Value")
ggsave(paste0(out_dir, "gene_list.gene_expression_violin_plot.p30.png"), width=8, height=4)
ggsave(paste0(out_dir, "gene_list.gene_expression_violin_plot.p30.svg"), width=8, height=4)



