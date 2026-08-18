library(CASSIA)
library(stringr)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/Harris-AM-6617/"

out_dir <- paste0(root_dir, "results/cassia_annotate/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))

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
  max_workers = 4                              # Number of parallel workers
)

# this cost 19 cents

# score the results
# 
# runCASSIA_score_batch(
#   input_file = paste0(out_dir, "results_gemini_lc_summary"),  # JSON auto-detected
#   output_file = paste0(out_dir, "results_gemini_lc_annotation_scored.csv"),
#   model = "anthropic/claude-sonnet-4.6",
#   provider = "openrouter"
# )


# try out claude for annotation (the recommended model)

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_claude_sonnet"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "anthropic/claude-sonnet-4.6",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)

# this cost $2.50


# # score it!
# runCASSIA_score_batch(
#   input_file = paste0(out_dir, "results_claude_sonnet_summary"),  # JSON auto-detected
#   output_file = paste0(out_dir, "results_claude_sonnet_annotation_scored.csv"),
#   model = "anthropic/claude-sonnet-4.6",
#   provider = "openrouter"
# )

# this was $1.50

