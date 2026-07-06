library(CASSIA)
library(gtools)
library(stringr)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ashley_bolte/human_meninges/"

out_dir <- paste0(root_dir, "results/cassia_annotate_total/")
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
  tissue = "Meninges",                             # Tissue type
  species = "Human",                           # Species
  model = "google/gemini-3-flash-preview",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)

# this cost 11 cents

# try out claude for annotation (the recommended model)

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_claude_sonnet"),   # Output file name
  tissue = "Meninges",                             # Tissue type
  species = "Human",                           # Species
  model = "anthropic/claude-sonnet-4.6",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)

# this cost $1.35




