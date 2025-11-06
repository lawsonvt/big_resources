# following the code from here https://chuiqin.github.io/irGSEA/
# tutorial also on github page https://github.com/chuiqin/irGSEA

library(SeuratData)
library(Seurat)
library(irGSEA)

# get data for tutorial
InstallData("pbmc3k")

data("pbmc3k.final")
# Update old Seurat object to accommodate new features
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)

DimPlot(pbmc3k.final, reduction = "umap",
        group.by = "seurat_annotations",label = T) + NoLegend()

Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations

# RUN ALL THE THINGS
# against the Hallmark dataset
pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 4,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

# creates assays for each method
Seurat::Assays(pbmc3k.final)

result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "seurat_annotations",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

summary(result.dge$ssgsea)

# NOTE: RRA (Robust Rank Aggreggation) aggregate the results across the different methods
names(result.dge)

head(result.dge$RRA)

rra_results <- result.dge$RRA

# heatmap
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot

# bubble plot
irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)
irGSEA.bubble.plot

# upset plot
irGSEA.upset.plot <- irGSEA.upset(object = result.dge, 
                                  method = "RRA")
irGSEA.upset.plot

# stacked bar plot
irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                      method = c("AUCell", "UCell", "singscore",
                                                 "ssgsea", "JASMINE", "viper", "RRA"))
irGSEA.barplot.plot


# umap scatterplot with geneset enrichment
scatterplot <- irGSEA.density.scatterplot(object = pbmc3k.final,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
                                          reduction = "umap")
scatterplot
irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "UCell", 
                                    top = 50)
irGSEA.bubble.plot


# cairo issues ... fixed when I installed X11 for mac

# violin plot and box plot
halfvlnplot <- irGSEA.halfvlnplot(object = pbmc3k.final,
                                  method = "UCell",
                                  show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
halfvlnplot

# compare across the various methods
vlnplot <- irGSEA.vlnplot(object = pbmc3k.final,
                          method = c("AUCell", "UCell", "singscore", "ssgsea", 
                                     "JASMINE", "viper"),
                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
vlnplot

# ridge plot
ridgeplot <- irGSEA.ridgeplot(object = pbmc3k.final,
                              method = "UCell",
                              show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
ridgeplot

irGSEA.ridgeplot(object = pbmc3k.final,
                 method = "ssgsea",
                 show.geneset = "HALLMARK-APOPTOSIS")


densityheatmap <- irGSEA.densityheatmap(object = pbmc3k.final,
                                        method = "UCell",
                                        show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
densityheatmap

# calculate the hub gene of the geneset

hub.result <- irGSEA.hub(object = pbmc3k.final, assay = "RNA", slot = "data",
                         method = c("AUCell","UCell","singscore", "ssgsea",
                                    "JASMINE", "viper"),
                         show.geneset = c("HALLMARK-INFLAMMATORY-RESPONSE",
                                          "HALLMARK-APOPTOSIS"),
                         ncores = 4, type = "rank", maxRank = 2000, top = 5,
                         correlation.color = c("#0073c2","white","#efc000"),
                         method.color = NULL)

head(hub.result$hub_result)

hub.result$hub_plot$`HALLMARK-APOPTOSIS`

# use cluster profiler to download other genesets

# try doing a GO analysis
pbmc3k.final.go <- pbmc3k.final

Seurat::Assays(pbmc3k.final)
# drop all previous analyses
pbmc3k.final.go@assays$AUCell <- NULL
pbmc3k.final.go@assays$UCell <- NULL
pbmc3k.final.go@assays$singscore <- NULL
pbmc3k.final.go@assays$ssgsea <- NULL
pbmc3k.final.go@assays$JASMINE <- NULL
pbmc3k.final.go@assays$viper <- NULL

# new analysis, focusing on GO terms
pbmc3k.final.go <- irGSEA.score(object = pbmc3k.final.go, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 4,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", 
                             category = "C5",  
                             subcategory = "GO:BP", 
                             geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')



