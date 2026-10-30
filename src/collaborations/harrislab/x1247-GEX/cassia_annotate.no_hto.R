library(CASSIA)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/x1247-GEX/"

out_dir <- paste0(root_dir, "results/cassia_annotate.no_hto/")
dir.create(out_dir, showWarnings = F)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/seurat_qc_samples.with_hto/all_markers.RDS"))

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_gemini_lc"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  model = "google/gemini-3-flash-preview",       # Model to use
  provider = "openrouter",                     # API provider
  max_workers = 4                              # Number of parallel workers
)


