library(CASSIA)
library(gtools)
library(stringr)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/cassia_endothelial_subclustering.pct_diff/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/endothelial_subclustering/all_markers.RDS"))

# rename clusters
levels(all_markers$cluster) <- paste0("cluster", str_pad(levels(all_markers$cluster), width=2, side="left", pad="0"))

# annotate!

# start with the cheapest model

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_gemini_lc"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "google/gemini-3-flash-preview",       # Model to use
  provider = "openrouter",                     # API provider
  ranking_method = "pct_diff",
  max_workers = 4                              # Number of parallel workers
)

# this cost 8 cents


# try out claude for annotation (the recommended model)

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_claude_sonnet"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "anthropic/claude-sonnet-4.6",       # Model to use
  provider = "openrouter",                     # API provider
  ranking_method = "pct_diff",
  max_workers = 4                              # Number of parallel workers
)

# this cost 16 cents


# # score the results
# 
# runCASSIA_score_batch(
#    input_file = paste0(out_dir, "results_gemini_lc_summary"),  # JSON auto-detected
#    output_file = paste0(out_dir, "results_gemini_lc_annotation_scored.csv"),
#    model = "anthropic/claude-sonnet-4.6",
#    provider = "openrouter"
# )
# 
# runCASSIA_score_batch(
#    input_file = paste0(out_dir, "results_claude_sonnet_summary"),  # JSON auto-detected
#    output_file = paste0(out_dir, "results_claude_sonnet_annotation_scored.csv"),
#    model = "anthropic/claude-sonnet-4.6",
#    provider = "openrouter"
# )


