library(Seurat)
library(SeuratObject)
library(CellChat)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(stringr)
library(snakecase)
library(dplyr)

root_dir <- "~/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/immune_cell_split/cellchat_analysis.immune/")
dir.create(out_dir, showWarnings = F)

# Read in and clean up data ----------------------------------------------------

# read in integrated seurat
seu_obj <- LoadSeuratRds(paste0(root_dir,
                                "results/celltype_naming/all_samples.celltype_named.seurat.RDS"))

# set assay to RNA and join layers
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- JoinLayers(seu_obj)

seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)

# load in immune subclustering metadata
immune_metadata <- readRDS(paste0(root_dir, "results/immune_cell_split/immune_cell_subclustering/subset_immune.metadata.RDS"))
nonimmune_metadata <- readRDS(paste0(root_dir, "results/immune_cell_split/nonimmune_cell_subclustering/subset_nonimmune.metadata.RDS"))

# add in cell IDs
immune_metadata$cell_id <- rownames(immune_metadata)
nonimmune_metadata$cell_id <- rownames(nonimmune_metadata)

# drop cell type (we gonna rename)
immune_metadata$cell_type <- NULL
immune_metadata$immune_clusters <- NULL

nonimmune_metadata$cell_type <- NULL

immune_cluster_raw <- read.xlsx(paste0(root_dir, "immune cell classification.xlsx"))

# reformat
immune_cluster_celltypes <- lapply(immune_cluster_raw$Cell.type, function(cell) {
  
  clusters <- unlist(strsplit(immune_cluster_raw[immune_cluster_raw$Cell.type == cell,]$Cluster, ","))
  
  return(data.frame(celltype=cell,
                    cluster=clusters))

})
immune_cluster_celltypes <- bind_rows(immune_cluster_celltypes)

# remove
immune_cluster_celltypes <- immune_cluster_celltypes[!grepl("Remove",immune_cluster_celltypes$celltype),]

# merge into metadata
immune_metadata <- merge(immune_metadata,
                         immune_cluster_celltypes,
                         by.x="harmony_clusters",
                         by.y="cluster")


# non immune
opc_cells <- as.character(c(2,10,18,19,20))

opc_metadata <-  nonimmune_metadata[nonimmune_metadata$nonimmune_clusters %in% opc_cells,]

opc_metadata$nonimmune_clusters <- NULL
opc_metadata$celltype <- "OPC"

# concat metadata
total_metadata <- rbind(immune_metadata,
                        opc_metadata)

rownames(total_metadata) <- total_metadata$cell_id

# because cell chat wants it
total_metadata$samples <- total_metadata$sample

# subset and redo metadata
seu_obj <- subset(seu_obj, cells=total_metadata$cell_id)

# assign new metadata
seu_obj@meta.data <- total_metadata[colnames(seu_obj),]

# fix idents
seu_obj@meta.data$celltype <- factor(seu_obj@meta.data$celltype,
                                     levels=sort(unique(seu_obj@meta.data$celltype)))

Idents(seu_obj) <- "celltype"


# Begin the cell chat analysis ------------------------------------------------

# split up Seurat

seu_wt <- subset(seu_obj, subset= condition == "WT")
seu_ko <- subset(seu_obj, subset= condition == "KO")

cellChat_list <- list(WT=createCellChat(object = seu_wt, group.by = "ident", assay = "RNA"),
                      KO=createCellChat(object = seu_ko, group.by = "ident", assay = "RNA"))

# remove all seurat objects and clean memory
rm(seu_obj)
rm(seu_wt)
rm(seu_ko)
gc()

# set the DB to be mouse
CellChatDB <- CellChatDB.mouse

showDatabaseCategory(CellChatDB)

# subset the Db to remove non-protein signaling
CellChatDB.use <- subsetDB(CellChatDB)

# add it to cell chat object
cellChat_list$WT@DB <- CellChatDB.use
cellChat_list$KO@DB <- CellChatDB.use

# process the cell chat objects
cellChat_list <- lapply(cellChat_list, function(cellChat) {
  
  # subset the expression data of signaling genes for saving computation cost
  cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 2) # do parallel
  
  options(future.globals.maxSize = 3 * 1e9) # increase the max future size from 500MB to 1GB
  
  cellChat <- identifyOverExpressedGenes(cellChat)
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  # smooth data (replaces projectData function, to apply the results to a PPI)
  # https://github.com/jinworks/CellChat/issues/185
  
  cellChat <- smoothData(cellChat, adj=PPI.mouse)
  
  # Compute the communication probability and infer cellular communication network
  cellChat <- computeCommunProb(cellChat, type = "triMean")
  
  # filter out low cell count communication
  cellChat <- filterCommunication(cellChat, min.cells = 10)
  
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellChat)
  
  # Calculate the aggregated cell-cell communication network
  cellChat <- aggregateNet(cellChat)
  
  
  
})

saveRDS(cellChat_list, file=paste0(out_dir, "cellchat_list.RDS"))

# Individual plots -------------------------------------------------------------

# overall interaction plots

for (condition in names(cellChat_list)) {
  
  cellChat <- cellChat_list[[condition]]
  
  groupSize <- as.numeric(table(cellChat@idents))
  
  pdf(paste0(out_dir, condition, ".number_of_interactions.circle_plot.pdf"), width=8, height=7)
  netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                   title.name = "Number of interactions")
  dev.off()
  
  pdf(paste0(out_dir, condition, ".weight_of_interactions.circle_plot.pdf"), width=8, height=7)
  netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                   title.name = "Interaction weights/strength")
  dev.off()
  
}

# cell based plots

for (condition in names(cellChat_list)) {
  
  cellChat <- cellChat_list[[condition]]
  
  groupSize <- as.numeric(table(cellChat@idents))
  
  cell_dir <- paste0(out_dir, condition, "_cell_plots/")
  
  dir.create(cell_dir, showWarnings = F)
  
  mat <- cellChat@net$weight
  
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    
    cell <- rownames(mat)[i]
    pdf(paste0(cell_dir, to_snake_case(cell), ".", condition, ".weight_of_interactions.circle_plot.pdf"), width=8, height=7)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
  }
  
  
}

# pathway based plots

# Infer the cell-cell communication at a signaling pathway level

cellChat_list <- lapply(cellChat_list, function(cellChat) {
  cellChat <- computeCommunProbPathway(cellChat)
  netAnalysis_computeCentrality(cellChat, slot.name = "netP")
})

for (condition in names(cellChat_list)) {
  
  cellChat <- cellChat_list[[condition]]
  
  groupSize <- as.numeric(table(cellChat@idents))
  
  pathway_dir <- paste0(out_dir, condition, "_pathway_plots/")
  
  dir.create(pathway_dir, showWarnings = F)
  
  sig_pathways <- cellChat@netP$pathways
  
  for (pathway in sig_pathways) {
    
    pdf(paste0(pathway_dir, pathway, ".", condition, ".chord_plot.pdf"), width=8, height=7)
    netVisual_aggregate(cellChat, signaling = pathway, 
                        layout = "chord")
    dev.off()
    
    pdf(paste0(pathway_dir, pathway, ".", condition, ".heatmap_plot.pdf"), width=8, height=7)
    print(netVisual_heatmap(cellChat, signaling = pathway, 
                            color.heatmap = 'Reds'))
    dev.off()
    
    pdf(paste0(pathway_dir, pathway, ".", condition, ".signaling_heatmap_plot.pdf"), width=8, height=7)
    print(netAnalysis_signalingRole_network(cellChat, 
                                            signaling = pathway, 
                                            width = 10, height = 4.5, font.size = 10))
    dev.off()
    
  }
  
}


saveRDS(cellChat_list, file=paste0(out_dir, "cellchat_list.RDS"))

# Comparison analysis ----------------------------------------------------------




