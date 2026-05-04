library(Seurat)
library(ggplot2)
library(stringr)
library(gtools)
library(cowplot)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/gene_expression_plots.sex_time_comp/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/find_marker_clusters/partially_celltype_named.seurat.RDS"))

# read in DEGs
total_degs <- readRDS(paste0(root_dir,
                             "results/differential_expression_analysis.sex_time_comp/total_results.RDS"))

# clean up metadata
metadata <- seu_obj@meta.data

metadata$sex <- factor(gsub("[Pp0-9]+-", "", metadata$MULTI_ID))
metadata$age <- toupper(gsub("-(female|male)", "", metadata$MULTI_ID))
metadata$age <- factor(metadata$age, levels=mixedsort(unique(metadata$age)))

metadata$sample_id <- str_to_title(metadata$MULTI_ID)
metadata$sample_id <- factor(metadata$sample_id,
                             levels=mixedsort(unique(metadata$sample_id)))

# add to seurat object
seu_obj$sample_id <- metadata$sample_id
seu_obj$sex <- metadata$sex
seu_obj$age <- metadata$age

# the cell types
cells <- c("Astrocyte",
           "Microglia")

# make subsets for plotting
astro_sub <- subset(seu_obj, subset = cell_type_partial == "Astrocyte")
micro_sub <- subset(seu_obj, subset = cell_type_partial == "Microglia")

rm(seu_obj)
gc()

total_degs$P14$Astrocyte

example_genes <- c("Xist","Uty")

p1 <- VlnPlot(astro_sub,
        features=example_genes[1],
        group.by = "age",
        split.by = "sex") +
  labs(x=NULL)

p2 <- VlnPlot(astro_sub,
              features=example_genes[2],
              group.by = "age",
              split.by = "sex") +
  labs(x=NULL)

plot_grid(p1,p2, ncol=1)
ggsave(paste0(out_dir, "example_violins.astrocyte.png"), width=6, height=5)

# custom plot: Mertk in microglia

VlnPlot(micro_sub,
        features="Mertk",
        group.by = "age",
        split.by = "sex") +
  labs(x=NULL)
ggsave(paste0(out_dir, "mertk_violins.microglia.png"), width=6, height=4)

# trim down further just to P14
astro_sub <- subset(astro_sub, subset = age == "P14")
micro_sub <- subset(micro_sub, subset = age == "P14")



# plot gene lists
micro_genes <- "Dlgap2, Cntnap4, Dlg2, Sorbs2, Cntn5, Gabra2, Grin2a, Grm8, Cntn4, Faah, Gabrb2, Cntn3, Nlgn1, Nrxn3, Cntnap2, Gabrb1, Cacna1c, Grip1, Grm7, Dlgap1, Adgrl3, Dlgap3, Grm5, Gria1, Shank2, Gria4, Adra1a, Gria3, Nrxn1, Grm1, Cntn6, Scn1a, Gabrb3, Cacng3, Prkaa2, Aprt, Grin3a, Scn3a, Pdgfra, Gria2, Gad2, Rit2, Shank1, Cntn1, Grm3, Nrxn2"
micro_genes <- trimws(unlist(strsplit(micro_genes, ", ")))

astro_genes <- "Grm7, Dlgap2, Cntn4, Cntnap2, Dlg2, Grm1, Cntn3, Gria1, Grm8, Cacna1c, Cntnap4, Gria4, Cntn5, Cacng3, Gad1"
astro_genes <- trimws(unlist(strsplit(astro_genes, ", ")))

# loop through
for (gene in micro_genes) {
  
  VlnPlot(micro_sub,
          features=gene,
          group.by = "age",
          split.by = "sex") +
    labs(x=NULL)
  ggsave(paste0(out_dir, gene, "_violins.microglia.png"), width=4, height=4)
  
}

for (gene in astro_genes) {
  
  VlnPlot(astro_sub,
          features=gene,
          group.by = "age",
          split.by = "sex") +
    labs(x=NULL)
  ggsave(paste0(out_dir, gene, "_violins.astrocytes.png"), width=4, height=4)
  
}





