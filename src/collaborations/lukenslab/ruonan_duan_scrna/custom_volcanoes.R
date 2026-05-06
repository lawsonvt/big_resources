library(ggplot2)
library(ggrepel)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/custom_volcanoes/")
dir.create(out_dir, showWarnings = F)


results <- readRDS(paste0(root_dir, "results/diff_exp_analysis.pseudobulk.cell_type/ko_minus_wt.cell_types.de_results.RDS"))

# volcano 1

cell <- "OPC"

genelist <- c("H2-D1",
              "H2-K1")


subset <- results[[cell]]

subset$log_p <- -log10(subset$pvalue)

subset_sig <- subset[subset$padj < 0.05 &
                       abs(subset$log2FoldChange) > 0.5,]


logp_thresh <- min(subset_sig$log_p)

ggplot(subset,
             aes(x=log2FoldChange,
                 y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  geom_hline(yintercept = logp_thresh,
             color="red", linetype=2) +
  geom_vline(xintercept = 0.5,
             color="red", linetype=2) +
  geom_vline(xintercept = -0.5,
             color="red", linetype=2) +
#  geom_point(data=subset_sig,
#             color="red") +
  geom_text_repel(data=subset[subset$gene %in% genelist,],
                  aes(label=gene),
                  color="blue", size=2.5,
                  max.overlaps = 50) +
  theme_bw() +
  labs(x="Log2 Fold Change", y="-log10(P-Value)", title=cell,
       subtitle="KO - WT")
ggsave(paste0(out_dir, "opc.volcanoes.H2_genes.png"), width=5, height=4)


