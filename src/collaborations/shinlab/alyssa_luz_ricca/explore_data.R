library(openxlsx)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(ggiraph)
library(patchwork)

root_dir <- "~/Documents/projects/shinlab/alyssa_luz_ricca/"

out_dir <- paste0(root_dir, "results/explore_data/")
dir.create(out_dir, showWarnings = F, recursive = T)

data <- read.xlsx(paste0(root_dir, "BAF_JS_4593_Proteins.xlsx"))

# clean up data

# remove data where all 3 abundance estimates are NA
data <- data[!(is.na(data$`Abundance:.F3:.Sample,.WT.Crtl`) &
               is.na(data$`Abundance:.F2:.Sample,.XIRP2/PTPRQ`) &
               is.na(data$`Abundance:.F1:.Sample,.XIRP2`)),]

# read in gene lists
ext_data_file <- paste0(root_dir, "GSE60019_gene_lists.xlsx")

sheet_names <- getSheetNames(ext_data_file)

gene_lists <- lapply(sheet_names, function(x) {
  
  
  raw_genes <- read.xlsx(ext_data_file, sheet = x, colNames=F)
  
  raw_genes <- unlist(raw_genes)
  raw_genes <- raw_genes[!is.na(raw_genes)]
  
  data.frame(type=x,
             gene_symbol=raw_genes)
  
  
})
gene_lists <- bind_rows(gene_lists)

gene_lists$cell_category <- "hair cell"
gene_lists[gene_lists$type == "surrounding",]$cell_category <- "surrounding cell" 

# merge them up

data <- merge(data,
              gene_lists,
              by.x="Gene.Symbol",
              by.y="gene_symbol",
              all.x=T)


data[is.na(data$`Abundance.Ratio.Adj..P-Value:.(XIRP2)./.(WT.Crtl)`),]$`Abundance.Ratio:.(XIRP2)./.(WT.Crtl)`
data[is.na(data$`Abundance.Ratio.Adj..P-Value:.(XIRP2/PTPRQ)./.(WT.Crtl)`),]$`Abundance.Ratio:.(XIRP2/PTPRQ)./.(WT.Crtl)`

# volcano plots

data$xirp2_pval <- data$`Abundance.Ratio.Adj..P-Value:.(XIRP2)./.(WT.Crtl)`
data[is.na(data$xirp2_pval),]$xirp2_pval <- 1

data$xirp2_ptprq_pval <- data$`Abundance.Ratio.Adj..P-Value:.(XIRP2/PTPRQ)./.(WT.Crtl)`
data[is.na(data$xirp2_ptprq_pval),]$xirp2_ptprq_pval <- 1


data$logp_xirp2_wt <- -log10(data$xirp2_pval)
data$logp_xirp2_ptprq_wt <- -log10(data$xirp2_ptprq_pval)

data$log2fc_xirp2_wt <- log2(data$`Abundance.Ratio:.(XIRP2)./.(WT.Crtl)`)
data$log2fc_xirp2_ptprq_wt <- log2(data$`Abundance.Ratio:.(XIRP2/PTPRQ)./.(WT.Crtl)`)

# the median isn't zero, shift it
data$log2fc_xirp2_wt <- data$log2fc_xirp2_wt - median(data$log2fc_xirp2_wt, na.rm=T)
data$log2fc_xirp2_ptprq_wt <- data$log2fc_xirp2_ptprq_wt - median(data$log2fc_xirp2_ptprq_wt, na.rm=T)

# volcano plots
v1 <- ggplot(data,
             aes(x=log2fc_xirp2_wt,
                 y=logp_xirp2_wt)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color="red", linetype=2) +
  labs(title="XIRP2 / WT", x="Log2(Abundance Ratio)", y="-Log10(Adj P-Value)")

v2 <- ggplot(data,
             aes(x=log2fc_xirp2_ptprq_wt,
                 y=logp_xirp2_ptprq_wt)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color="red", linetype=2) +
  labs(title="(XIRP2/PTPRQ) / WT", x="Log2(Abundance Ratio)", y="-Log10(Adj P-Value)")

plot_grid(v1, v2)

ggsave(paste0(out_dir, "volcano_plots.png"), width=8, height=5)

# label significance
data$sig_type <- as.numeric(data$xirp2_pval < 0.05) +
  as.numeric(data$xirp2_ptprq_pval < 0.05) *2
data$sig_type <- factor(data$sig_type)
levels(data$sig_type) <- c("None","XIRP2", "XIRP2/PTPRQ", "Both")

data <- data[order(data$sig_type),]

#View(data[is.na(data$sig_type),])


ggplot(data,
       aes(x=log2fc_xirp2_wt,
           y=log2fc_xirp2_ptprq_wt,
           color=sig_type)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey") +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("grey",
                              "blue",
                              "red",
                              "purple")) +
  theme_bw() +
  labs(x="Log2FC XIRP2 / WT", y="Log2FC (XIRP2/PTPRQ) / WT",
       color="Adj P < 0.05")
ggsave(paste0(out_dir, "foldchange_scatter_plot.total_proteins.png"), width=7, height=5)

ggplot(data,
       aes(x=log2fc_xirp2_wt,
           y=log2fc_xirp2_ptprq_wt,
           color=sig_type)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey") +
  facet_wrap(~ cell_category, ncol=3) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("grey",
                              "blue",
                              "red",
                              "purple")) +
  theme_bw() +
  labs(x="Log2FC XIRP2 / WT", y="Log2FC (XIRP2/PTPRQ) / WT",
       color="Adj P < 0.05")
ggsave(paste0(out_dir, "foldchange_scatter_plot.celltype_proteins.png"), width=12, height=5)


# drop temp pvals
data$xirp2_pval <- NULL
data$xirp2_ptprq_pval <- NULL

# output to an excel file

data <- data[order(data$sig_type, decreasing=T),]

export <- list(total=data,
               no_surrounding=data[is.na(data$cell_category) | data$cell_category != "surrounding cell",],
               hair_cell=data[!is.na(data$cell_category) & data$cell_category == "hair cell",])

write.xlsx(export, paste0(out_dir, "categorized_data.xlsx"))

# make an interactive plot?

data <- data[order(data$sig_type),]

p <- ggplot(data,
       aes(x=log2fc_xirp2_wt,
           y=log2fc_xirp2_ptprq_wt,
           color=sig_type,
           tooltip = Gene.Symbol)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey") +
  geom_point_interactive(alpha=0.5) +
  scale_color_manual(values=c("grey",
                              "blue",
                              "red",
                              "purple")) +
  theme_bw() +
  labs(x="Log2FC XIRP2 / WT", y="Log2FC (XIRP2/PTPRQ) / WT",
       color="Adj P < 0.05")

int_p <- girafe(p, width_svg = 7, height_svg = 5)
htmltools::save_html(int_p, paste0(out_dir, "foldchange_scatter_plot.total_proteins.html"))

p <- ggplot(data,
       aes(x=log2fc_xirp2_wt,
           y=log2fc_xirp2_ptprq_wt,
           color=sig_type,
           tooltip=Gene.Symbol)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey") +
  facet_wrap(~ cell_category, ncol=3) +
  geom_point_interactive(alpha=0.5) +
  scale_color_manual(values=c("grey",
                              "blue",
                              "red",
                              "purple")) +
  theme_bw() +
  labs(x="Log2FC XIRP2 / WT", y="Log2FC (XIRP2/PTPRQ) / WT",
       color="Adj P < 0.05")

int_p <- girafe(p, width_svg = 12, height_svg = 5)


htmltools::save_html(int_p,
                     paste0(out_dir, "foldchange_scatter_plot.celltype_proteins.html"))

# interactive volcano plots

v1 <- ggplot(data,
             aes(x=log2fc_xirp2_wt,
                 y=logp_xirp2_wt,
                 tooltip=Gene.Symbol,
                 data_id=Gene.Symbol)) +
  geom_point_interactive(alpha=0.5) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color="red", linetype=2) +
  labs(title="XIRP2 / WT", x="Log2(Abundance Ratio)", y="-Log10(Adj P-Value)")

v2 <- ggplot(data,
             aes(x=log2fc_xirp2_ptprq_wt,
                 y=logp_xirp2_ptprq_wt,
                 tooltip=Gene.Symbol,
                 data_id=Gene.Symbol)) +
  geom_point_interactive(alpha=0.5) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color="red", linetype=2) +
  labs(title="(XIRP2/PTPRQ) / WT", x="Log2(Abundance Ratio)", y="-Log10(Adj P-Value)")

combined_plot <- v1 + v2 + plot_layout(ncol=2)

int_p <- girafe(combined_plot, width_svg=8, height_svg=5)

# Set options for the interactive plot
int_p <- girafe_options(
  int_p,
  opts_hover(css = "fill:cyan;stroke:black;cursor:pointer;"),
  opts_selection(type = "single", css = "fill:red;stroke:black;")
)

htmltools::save_html(int_p,
                     paste0(out_dir, "volcano_plots.html"))



