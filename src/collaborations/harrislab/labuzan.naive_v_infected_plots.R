library(ggplot2)
library(ggrepel)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)

out_dir <- "~/Documents/projects/harrislab/naive_v_infected_brain_mouse/R/plots/"

dir.create(out_dir, showWarnings = F)

degs <- read.csv("~/Documents/projects/harrislab/naive_v_infected_brain_mouse/R/gene_level_results.csv")

# filter out NAs
degs <- degs[!is.na(degs$padj),]

# order them
degs <- degs[order(degs$padj),]

# filter down to expressed genes
degs <- degs[degs$Naive_FPKM_Mean > 1 |
                 degs$Infected_FPKM_mean > 1,]


nrow(degs[degs$padj < 0.01,])

# add in data for plotting
degs$log_p <- - log10(degs$pvalue)

# replace Inf with max
max_logp <- max(degs[!is.infinite(degs$log_p),]$log_p)

degs[is.infinite(degs$log_p),]$log_p <- max_logp

# get GO terms
gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")

ag_pathway <- gobp_gene_sets[gobp_gene_sets$gs_exact_source == "GO:0019882",]

sig_degs <- degs[degs$padj < 0.01 &
                   abs(degs$log2FoldChange) > 0.5,]

sig_ag <- degs[degs$padj < 0.01 &
                 degs$ensembl_gene_id %in% ag_pathway$ensembl_gene &
                 abs(degs$log2FoldChange) > 0.5,]

logpthresh <- min(sig_degs$log_p)

ggplot(degs,
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.1) +
  geom_point(data=sig_degs, color="red", alpha=0.2) +
  geom_text_repel(data=sig_ag, color="black", aes(label=external_gene_name), 
                  max.overlaps = 50, size = 3) +
  theme_bw() +
  geom_vline(xintercept = 0.5, linetype=2) +
  geom_vline(xintercept = -0.5, linetype=2) +
  geom_hline(yintercept = logpthresh, linetype=2) +
  labs(x="Log2 Fold Change", y="-Log10 P")
ggsave(paste0(out_dir, "volcano_plot.antigen_pres_annotate.png"), width=7, height=6)

# heatmap ----------------------------------------------------------------------

vst <- read.csv("~/Documents/projects/harrislab/naive_v_infected_brain_mouse/R/gene_level_vst_counts.csv")

vst_sig_ag <- vst[vst$ensembl_gene_id %in% sig_ag$ensembl_gene_id,]
rownames(vst_sig_ag) <- vst_sig_ag$external_gene_name

# make into matrix
vst_sig_ag$external_gene_name <- NULL
vst_sig_ag$ensembl_gene_id <- NULL

vst_sig_ag <- as.matrix(vst_sig_ag)

# scale it
vst_sig_ag_s <- t(scale(t(vst_sig_ag)))

pdf(paste0(out_dir, "heatmap.antigen_pres.pdf"), width=6, height=9)
print(Heatmap(vst_sig_ag_s,
        name = "VST\nZ-Score",
        col=colorRamp2(c(-1.5, 0, 1.5), c("blue","white","red")),
        row_names_gp = gpar(fontsize = 9)))
dev.off()

