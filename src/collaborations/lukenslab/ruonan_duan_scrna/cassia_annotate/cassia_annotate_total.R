library(CASSIA)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/lukenslab/ruonan_duan/"
#root_dir <- "~/projects/lukenslab/ruonan_duan/"

out_dir <- paste0(root_dir, "results/cassia_annotate/cassia_annotate_total/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/cassia_annotate/find_all_markers.total/all_markers.RDS"))

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

# this cost 25 cents

# score the results

# (tried using open ai GPT 5.4, but kept getting http errors and it wouldnt fully score the cluster names)

runCASSIA_score_batch(
  input_file = paste0(out_dir, "results_gemini_lc_summary"),  # JSON auto-detected
  output_file = paste0(out_dir, "results_gemini_lc_annotation_scored.csv"),
  model = "anthropic/claude-sonnet-4.6",
  provider = "openrouter"
)

# this cost around $1


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

# this cost over $2


# score it!
runCASSIA_score_batch(
  input_file = paste0(out_dir, "results_claude_sonnet_summary"),  # JSON auto-detected
  output_file = paste0(out_dir, "results_claude_sonnet_annotation_scored.csv"),
  model = "anthropic/claude-sonnet-4.6",
  provider = "openrouter"
)

# this was $1.50

# OLLAMA

runCASSIA_batch(
  marker = all_markers,                # Marker data from FindAllMarkers
  output_name = paste0(out_dir, "results_ollama_gpt"),   # Output file name
  tissue = "Brain",                             # Tissue type
  species = "Mouse",                           # Species
  provider = "http://localhost:11434/v1",
  model = "gpt-oss:20b",             # API provider
  max_workers = 4                              # Number of parallel workers
)

