library(CASSIA)
library(openxlsx)
library(gtools)
library(Seurat)
library(SeuratObject)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"


out_dir <- paste0(root_dir, "results/cassia_results_analysis/")
dir.create(out_dir, showWarnings = F, recursive = T)


# read in CASSIA results
gemini_results <- read.csv(paste0(root_dir, "results/cassia_annotate/results_gemini_lc_summary.csv"))
claude_results <- read.csv(paste0(root_dir, "results/cassia_annotate/results_claude_sonnet_summary.csv"))
olgpt_results <- read.csv(paste0(root_dir, "results/cassia_annotate/results_ollama_gpt_summary.csv"))

# create smaller versions
gemini_simple <- gemini_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                                   "Predicted.Detailed.Cell.Type",
                                   "Possible.Mixed.Cell.Types")]
colnames(gemini_simple) <- c("harmony_clusters",
                             "gemini_prediction",
                             "gemini_detailed_prediction",
                             "gemini_possible_mix")

claude_simple <- claude_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                                   "Predicted.Detailed.Cell.Type",
                                   "Possible.Mixed.Cell.Types")]

colnames(claude_simple) <- c("harmony_clusters",
                             "claude_prediction",
                             "claude_detailed_prediction",
                             "claude_possible_mix")

olgpt_simple <- olgpt_results[,c("Cluster.ID","Predicted.General.Cell.Type",
                          "Predicted.Detailed.Cell.Type",
                          "Possible.Mixed.Cell.Types")]


colnames(olgpt_simple) <- c("harmony_clusters",
                             "ollama_prediction",
                             "ollama_detailed_prediction",
                             "ollama_possible_mix")




# merge em up for a comparison


cell_comp <- merge(claude_simple,
                   gemini_simple,
                   by="harmony_clusters")

cell_comp <- merge(cell_comp,
                   olgpt_simple,
                   by="harmony_clusters",
                   all.x=T)

cell_comp <- cell_comp[mixedorder(cell_comp$harmony_clusters),]


# output the comp

wb <- createWorkbook()
addWorksheet(wb, "comparison")
writeData(wb, "comparison", cell_comp)

n_rows <- nrow(cell_comp) + 1  # +1 for header row

setColWidths(wb, "comparison", cols=1, width=20)
setColWidths(wb, "comparison", cols=2:ncol(cell_comp), width=40)

setRowHeights(wb, "comparison", rows = 2:n_rows, heights = 80)

wrap_style <- createStyle(wrapText = TRUE)
addStyle(wb, "comparison",
         style = wrap_style,
         rows = 2:n_rows,
         cols = 2:ncol(cell_comp),       # whichever columns should wrap
         gridExpand = TRUE)    # applies style to every row/col combination

saveWorkbook(wb, paste0(out_dir, "cassia_results_comparison.xlsx"), overwrite = T)

saveRDS(cell_comp, paste0(out_dir, "cassia_results_comparison.RDS"))

# move report htmls into this folder

html_reports <- list.files(paste0(root_dir, "results/cassia_annotate/"),
                           "_report.html", full.names = T)
html_reports <- html_reports[!grepl("scored", html_reports)]

file.copy(from=html_reports, to=out_dir)

# look at UMAP

# read in integrated seurat
int_seu <- LoadSeuratRds(paste0(root_dir,
                                "results/integrate_seurat_samples/all_samples.integrated_seurat.RDS"))


DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters", label=T)


# determine the endothelial cells

endo_cells <- cell_comp[grepl("[Ee]ndothelial", cell_comp$gemini_prediction),]

endo_clusters <- gsub("cluster", "", endo_cells$harmony_clusters)

# subset the seurat
endo_seu <- subset(int_seu, subset = harmony_clusters %in% endo_clusters)

# clean it up
umap_coords <- Embeddings(endo_seu, reduction = "umap.harmony")
umap_coords <- as.data.frame(umap_coords)

outliers <- umap_coords[umap_coords$umapharmony_1 < 0 |
                          umap_coords$umapharmony_2 < -5,]

endo_seu$cell_id <- colnames(endo_seu)

endo_seu <- subset(endo_seu, subset = !cell_id %in% rownames(outliers))

DimPlot(endo_seu, reduction="umap.harmony", group.by= "harmony_clusters", label=T, raster = F, label.box=T)



