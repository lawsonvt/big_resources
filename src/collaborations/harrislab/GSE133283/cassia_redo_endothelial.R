library(CASSIA)
library(gtools)
library(stringr)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/cassia_redo_endothelial/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))

# rename clusters
levels(all_markers$cluster) <- paste0("cluster", str_pad(levels(all_markers$cluster), width=2, side="left", pad="0"))

# get previous cassia results
cassia_results <- readRDS(paste0(root_dir, "results/cassia_results_analysis/cassia_results_comparison.RDS"))

# subset the endothelial cells
endo_cells <- cassia_results[grepl("[Ee]ndothelial", cassia_results$gemini_prediction),]

endo_cells$harmony_clusters <- paste0("cluster",str_pad(gsub("cluster", "", endo_cells$harmony_clusters), width=2, side="left", pad="0"))

# subset the markers
endo_markers <- all_markers[all_markers$cluster %in% endo_cells$harmony_clusters,]

# redo factor
endo_markers$cluster <- factor(as.character(endo_markers$cluster),
                               levels=mixedsort(unique(as.character(endo_markers$cluster))))

# annotate!

# start with the cheapest model

runCASSIA_batch(
  marker = endo_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_gemini_lc"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "google/gemini-3-flash-preview",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)

# this cost 4 cents


# try out claude for annotation (the recommended model)

runCASSIA_batch(
  marker = endo_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_claude_sonnet"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "anthropic/claude-sonnet-4.6",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)

# this cost 49 cents


# score the results

runCASSIA_score_batch(
   input_file = paste0(out_dir, "results_gemini_lc_summary"),  # JSON auto-detected
   output_file = paste0(out_dir, "results_gemini_lc_annotation_scored.csv"),
   model = "anthropic/claude-sonnet-4.6",
   provider = "openrouter"
)

runCASSIA_score_batch(
   input_file = paste0(out_dir, "results_claude_sonnet_summary"),  # JSON auto-detected
   output_file = paste0(out_dir, "results_claude_sonnet_annotation_scored.csv"),
   model = "anthropic/claude-sonnet-4.6",
   provider = "openrouter"
)


