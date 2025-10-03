# following the code from here 
# https://www.borch.dev/uploads/screpertoire/articles/running_escape
# or not, unclear where the starting point is
# this is a better one 
# https://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/escape.html

library(escape) # install via BioConductor
library(scran)
library(Seurat)
library(SeuratObject)
library(RColorBrewer)
library(ggplot2)
library(BiocParallel)

# get example dataset from Seurat
pbmc_small <- get("pbmc_small")

# first step is to get gene sets
# multiple examples on how to do so

# this pulls from msigdb
# go here to get info on library names https://www.gsea-msigdb.org/gsea/msigdb
gs_hallmark <- getGeneSets(library="H", 
                           species="Homo sapiens") # also supports mouse

gs_c2 <- getGeneSets(library="C2",
                     subcategory = c("CP:BioCarta",
                                     "CP:KEGG_LEGACY",
                                     "CP:Reactome",
                                     "CP:WikiPathways"),
                     species="Homo sapiens")

# gs_go <- getGeneSets(library="C5",
#                      subcategory = c("GO:BP",
#                                      "GO:MF",
#                                      "GO:CC"),
#                      species = "Homo sapiens")

# built in gene sets
data("escape.gene.sets", package="escape")
gene_sets <- escape.gene.sets

# build your own gene sets!
custom_gene_sets <- list(Bcells = c("MS4A1","CD79B","CD79A","IGH1","IGH2"),
                  Myeloid = c("SPI1","FCER1G","CSF1R"),
                  Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))

# two main ways to run

# 1. escape.matrix generates a matrix of values: geneset scores per cell

enrichment_scores <- escape.matrix(pbmc_small, 
                                   method = "ssGSEA",
                                   gene.sets = gs_hallmark, 
                                   groups = 1000, 
                                   min.size = 5)
enrichment_scores[1:5,1:5]
dim(enrichment_scores)

# output will automatically remove genesets that do not have enough gene coverage

# example plot for vignette ... not sure what this is supposed to show but
# it looks the same
# ggplot(data = as.data.frame(enrichment_scores), 
#        mapping = aes(enrichment_scores[,1], enrichment_scores[,2])) + 
#   geom_point() + 
#   theme_classic() + 
#   theme(axis.title = element_blank())

# run in parallel to make more efficient
enrichment_scores <- escape.matrix(pbmc_small, 
                                   method = "ssGSEA",
                                   gene.sets = gs_hallmark, 
                                   groups = 1000, 
                                   min.size = 5,
                                   BPPARAM = SnowParam(workers = 2))

# 2. runEscape performs the same calculation and attaches it to the Seurat object
# (it can also be run in parallel)

pbmc_small <- runEscape(pbmc_small,
                        method = "ssGSEA",
                        gene.sets = gs_hallmark, 
                        groups = 1000, 
                        min.size = 5,
                        new.assay.name = "hallmark_ssGSEA")

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

DefaultAssay(pbmc_small) <- "hallmark_ssGSEA"
FeaturePlot(pbmc_small, "HALLMARK-APOPTOSIS") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank())

pbmc_small <- runEscape(pbmc_small,
                        method = "UCell",
                        gene.sets = gs_hallmark, 
                        groups = 1000, 
                        min.size = 5,
                        new.assay.name = "hallmark_UCell")

DefaultAssay(pbmc_small) <- "hallmark_UCell"
FeaturePlot(pbmc_small, "HALLMARK-APOPTOSIS") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank())

DefaultAssay(pbmc_small) <- "RNA"


# the following code doesnt run due to
# "Fewer than 20% of the genes in the gene sets are included in the rankings.
# Check wether the gene IDs in the 'rankings' and 'geneSets' match.

# pbmc_small <- runEscape(pbmc_small,
#                         method = "AUCell",
#                         gene.sets = gs_hallmark, 
#                         groups = 1000, 
#                         min.size = 5,
#                         new.assay.name = "hallmark_AUCell")

# but it works with a different geneset list
pbmc_small <- runEscape(pbmc_small,
                        method = "AUCell",
                        gene.sets = gene_sets, 
                        groups = 1000, 
                        min.size = 5,
                        new.assay.name = "escape_gs_AUCell")

pbmc_small <- runEscape(pbmc_small,
                        method = "GSVA",
                        gene.sets = gs_hallmark, 
                        groups = 1000, 
                        min.size = 5,
                        new.assay.name = "hallmark_GSVA")

# ok, the genes have been transformed to genesets, now what?

# we normalize!
pbmc_small <- performNormalization(sc.data=pbmc_small,
                                   assay="hallmark_ssGSEA",
                                   gene.sets = gs_hallmark)
# creates a new assay called 'hallmark_ssGSEA_normalized'

pbmc_small@assays$hallmark_ssGSEA_normalized$data[1:5,1:5]

# can also normalize based on the pre-calculated scaling factor
pbmc_small <- performNormalization(sc.data=pbmc_small,
                                   assay="hallmark_ssGSEA",
                                   gene.sets = gs_hallmark,
                                   scale.factor = pbmc_small$nFeature_RNA)

pbmc_small@assays$hallmark_ssGSEA_normalized$data[1:5,1:5]

# in addition, we can make the values positive for better log2FC
pbmc_small <- performNormalization(sc.data=pbmc_small,
                                   assay="hallmark_ssGSEA",
                                   gene.sets = gs_hallmark,
                                   make.positive = T)

pbmc_small@assays$hallmark_ssGSEA_normalized$data[1:5,1:5]

# Visualizations
heatmapEnrichment(pbmc_small, 
                  group.by = "RNA_snn_res.1", # metadata field
                  gene.set.use = "all",
                  assay = "hallmark_ssGSEA")

# add clustering (but no dendograms)
heatmapEnrichment(pbmc_small, 
                  group.by = "RNA_snn_res.1", # metadata field
                  gene.set.use = "all",
                  assay = "hallmark_ssGSEA",
                  cluster.rows = T,
                  cluster.columns = T)

# can use other palettes
hcl.pals()
heatmapEnrichment(pbmc_small,
                  assay="hallmark_ssGSEA",
                  palette="Spectral")

# alternatively, since the heatmap is a ggplot object, can alter with ggplot functions
heatmapEnrichment(pbmc_small,
                  assay="hallmark_ssGSEA") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 

# geyser plots
geyserEnrichment(pbmc_small, 
                 assay = "hallmark_ssGSEA",
                 gene.set = "HALLMARK-INTERFERON-GAMMA-RESPONSE")

# order by the mean value
geyserEnrichment(pbmc_small, 
                 assay = "hallmark_ssGSEA",
                 gene.set = "HALLMARK-INTERFERON-GAMMA-RESPONSE",
                 order.by = "mean")
# use facets
geyserEnrichment(pbmc_small, 
                 assay = "hallmark_ssGSEA",
                 gene.set = "HALLMARK-INTERFERON-GAMMA-RESPONSE", 
                 facet.by = "groups")

# color by y instead of x
geyserEnrichment(pbmc_small, 
                 assay = "hallmark_ssGSEA",
                 gene.set = "HALLMARK-INTERFERON-GAMMA-RESPONSE",
                 color.by = "HALLMARK-INTERFERON-GAMMA-RESPONSE")

# ridge plots!
ridgeEnrichment(pbmc_small, 
                assay = "hallmark_ssGSEA",
                gene.set = "HALLMARK-IL2-STAT5-SIGNALING")

# scale values and add locations of cells along x-axis (rug)
ridgeEnrichment(pbmc_small, 
                assay = "hallmark_ssGSEA",
                gene.set = "HALLMARK-IL2-STAT5-SIGNALING",
                add.rug = TRUE,
                scale = TRUE)




