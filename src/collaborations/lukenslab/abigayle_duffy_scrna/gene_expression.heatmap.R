library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(gtools)
library(stringr)
library(ggplot2)

root_dir <- "~/Documents/projects/lukenslab/abigayle_duffy/"

out_dir <- paste0(root_dir, "results/gene_expression.heatmap/")
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/find_marker_clusters/partially_celltype_named.seurat.RDS"))

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

# filter down seurat object and aggregate

seu_subset <- subset(seu_obj, subset = cell_type_partial %in% cells &
                        age == "P14")

# aggregate values

seu_pseudo <- AggregateExpression(
  seu_subset,
  assays = "RNA",
  slot = "counts",
  group.by = c("cell_type_partial","sex")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# metadata
metadata <- data.frame(cell=sapply(colnames(pseudobulk_matrix), function(x) {unlist(strsplit(x, "_"))[1]}),
                      sex=sapply(colnames(pseudobulk_matrix), function(x) {unlist(strsplit(x, "_"))[2]}))
rownames(metadata) <- colnames(pseudobulk_matrix)


# make DESeq object

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_matrix,
  colData = metadata,
  design = ~ cell + sex
)

# pull out vsd
vsd <- vst(dds, blind=F)

counts <- assay(vsd)


# plot gene lists
micro_genes <- "Dlgap2, Cntnap4, Dlg2, Sorbs2, Cntn5, Gabra2, Grin2a, Grm8, Cntn4, Faah, Gabrb2, Cntn3, Nlgn1, Nrxn3, Cntnap2, Gabrb1, Cacna1c, Grip1, Grm7, Dlgap1, Adgrl3, Dlgap3, Grm5, Gria1, Shank2, Gria4, Adra1a, Gria3, Nrxn1, Grm1, Cntn6, Scn1a, Gabrb3, Cacng3, Prkaa2, Aprt, Grin3a, Scn3a, Pdgfra, Gria2, Gad2, Rit2, Shank1, Cntn1, Grm3, Nrxn2"
micro_genes <- trimws(unlist(strsplit(micro_genes, ", ")))

astro_genes <- "Grm7, Dlgap2, Cntn4, Cntnap2, Dlg2, Grm1, Cntn3, Gria1, Grm8, Cacna1c, Cntnap4, Gria4, Cntn5, Cacng3, Gad1"
astro_genes <- trimws(unlist(strsplit(astro_genes, ", ")))

# make microglia heatmap
subset_mat <- counts[micro_genes,c("Microglia_female","Microglia_male")]
colnames(subset_mat) <- c("female", "male")

subset_mat <- t(subset_mat)

pdf(paste0(out_dir, "microglia_p14.sex_comp_heatmap.pdf"),
    width=10, height=4)
print(Heatmap(subset_mat,
              name = "VST Expr"))
dev.off()

# make astrocyte heatmap

subset_mat <- counts[astro_genes,c("Astrocyte_female","Astrocyte_male" )]
colnames(subset_mat) <- c("female", "male")

subset_mat <- t(subset_mat)

pdf(paste0(out_dir, "astrocyte_p14.sex_comp_heatmap.pdf"),
    width=8, height=4)
print(Heatmap(subset_mat,
              name = "VST Expr"))
dev.off()

# alternative idea: subset dotplots
seu_subset_micro <- subset(seu_subset, subset = cell_type_partial == "Microglia")
seu_subset_astro <- subset(seu_subset, subset = cell_type_partial == "Astrocyte")

Idents(seu_subset_micro) <- "sex"
Idents(seu_subset_astro) <- "sex"

DotPlot(seu_subset_micro,
        features = micro_genes,
        scale=F) + RotatedAxis() + labs(x=NULL, y=NULL, title="Microglia P14 Cells")
ggsave(paste0(out_dir, "microglia_p14.sex_comp_dotplot.png"), width=12, height=4, bg="white")

DotPlot(seu_subset_astro,
        features = astro_genes,
        scale=F) + RotatedAxis() + labs(x=NULL, y=NULL, title="Astrocyte P14 Cells")
ggsave(paste0(out_dir, "astrocyte_p14.sex_comp_dotplot.png"), width=10, height=4, bg="white")


# new genelist sets, for the P7 time point

p7_subset <- subset(seu_obj, subset = cell_type_partial %in% c("Astrocyte",
                                                             "Microglia") &
                       age == "P7")
p7_astro <- subset(p7_subset, subset = cell_type_partial == "Astrocyte")
p7_micro <- subset(p7_subset, subset = cell_type_partial == "Microglia")


astro_genes2 <- "Uty, Kdm5d, Hmgb2, Taf7, H1f0, H3f3b, Ing5, Zdbf2, Nsmce1, Msh6, Tspyl2, Lmnb1, Tspyl5, Apobec3, Suv39h1, Pwp1, Rbl1, Cxxc1, Thap7, Hells, Fbxo30, Mbd3, Pax6"
astro_genes2 <- trimws(unlist(strsplit(astro_genes2, ", ")))

micro_genes2 <- "Tonsl, Dtl, Gins1, Orc1, Gins2, Cdt1, Donson, E2f8, Rrm1, Mcm6, Pole, Chtf18, Atad5, Wdhd1, Tk1, Blm, Rfc3, Cdk1, Mcm10, Pola1, Gmnn, Mcm2, Pole2, Rfc2, Brca2, Pola2, Mcm4, Lig1, Fgfr1, Rmi2, Topbp1, Dbf4, Pold1, Mms22l, Gins4, Rpa2, Wee1, Exo1, Timeless, Rad51, Mcm5, Dna2, Ing5, Chaf1b, Actl6a, Fancm, Pcna, Prim2, Helb, Rnaseh2a, Rtel1, Orc6, Rev1, Ccne1, Fbxo5, Jade3, Nasp, Pclaf, Recql5, Rbms1, Samhd1, Zfp830, Prim1"
micro_genes2 <- trimws(unlist(strsplit(micro_genes2, ", ")))

# aggregate expression for the heatmap

# aggregate values

seu_pseudo <- AggregateExpression(
  p7_subset,
  assays = "RNA",
  slot = "counts",
  group.by = c("cell_type_partial","sex")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# metadata
metadata <- data.frame(cell=sapply(colnames(pseudobulk_matrix), function(x) {unlist(strsplit(x, "_"))[1]}),
                       sex=sapply(colnames(pseudobulk_matrix), function(x) {unlist(strsplit(x, "_"))[2]}))
rownames(metadata) <- colnames(pseudobulk_matrix)


# make DESeq object

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_matrix,
  colData = metadata,
  design = ~ cell + sex
)

# pull out vsd
vsd <- vst(dds, blind=F)

counts <- assay(vsd)

# make microglia heatmap
subset_mat <- counts[micro_genes2,c("Microglia_female","Microglia_male")]
colnames(subset_mat) <- c("female", "male")

subset_mat <- t(subset_mat)

pdf(paste0(out_dir, "microglia_p7.sex_comp_heatmap.pdf"),
    width=12, height=4)
print(Heatmap(subset_mat,
              name = "VST Expr"))
dev.off()

# make astrocyte heatmap

subset_mat <- counts[astro_genes2,c("Astrocyte_female","Astrocyte_male" )]
colnames(subset_mat) <- c("female", "male")

subset_mat <- t(subset_mat)

pdf(paste0(out_dir, "astrocyte_p7.sex_comp_heatmap.pdf"),
    width=8, height=4)
print(Heatmap(subset_mat,
              name = "VST Expr"))
dev.off()


# dot plots
Idents(p7_micro) <- "sex"
Idents(p7_astro) <- "sex"

DotPlot(p7_micro,
        features = micro_genes2,
        scale=F) + RotatedAxis() + labs(x=NULL, y=NULL, title="Microglia P7 Cells")
ggsave(paste0(out_dir, "microglia_p7.sex_comp_dotplot.png"), width=14, height=4, bg="white")

DotPlot(p7_astro,
        features = astro_genes2,
        scale=F) + RotatedAxis() + labs(x=NULL, y=NULL, title="Astrocyte P7 Cells")
ggsave(paste0(out_dir, "astrocyte_p7.sex_comp_dotplot.png"), width=10, height=4, bg="white")





