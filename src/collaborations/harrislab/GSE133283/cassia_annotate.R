library(CASSIA)

validate_api_keys()

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/GSE133283/"

out_dir <- paste0(root_dir, "results/cassia_annotate/")
dir.create(out_dir, showWarnings = F, recursive = T)

# load in markers
all_markers <- readRDS(paste0(root_dir, "results/find_all_markers/all_markers.RDS"))

# rename clusters
all_markers$cluster <- paste0("cluster", all_markers$cluster)

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


# still having issues despite using 3 different models ...

# pull in results
ollama_gpt_results <- read.csv(paste0(out_dir, "results_ollama_gpt_summary.csv"))

total_clusters <- unique(all_markers$cluster)

missing_clusters <- total_clusters[!total_clusters %in% ollama_gpt_results$Cluster.ID]

all_markers <- all_markers[order(all_markers$avg_log2FC, decreasing=T),]

top_markers_missing <- all_markers[all_markers$cluster == missing_clusters[1] &
                                     all_markers$pct.1 > 0.1 &
                                     all_markers$avg_log2FC > 0.25 &
                                     all_markers$p_val_adj < 0.05,][1:50,]


result <- runCASSIA(
  marker_list = top_markers_missing$gene,
  provider = "http://localhost:11434/v1",
  model = "gpt-oss:20b",
  tissue = "Brain",                 
  species = "Mouse" 
)

# make a loop for the missing clusters

missing_results <- lapply(missing_clusters, function(cluster) {
  
  print(cluster)
  
  top_markers <- all_markers[all_markers$cluster == cluster &
                               all_markers$pct.1 > 0.1 &
                               all_markers$avg_log2FC > 0.25 &
                               all_markers$p_val_adj < 0.05,][1:50,]
  
  result <- runCASSIA(
    marker_list = top_markers$gene,
    provider = "http://localhost:11434/v1",
    model = "gpt-oss:20b",
    tissue = "Brain",                 
    species = "Mouse" 
  )
  
  return(result)
})

lapply(missing_results, function(result) {result$main_cell_type})

# try out a different model

missing_results <- lapply(missing_clusters, function(cluster) {
  
  print(cluster)
  
  top_markers <- all_markers[all_markers$cluster == cluster &
                               all_markers$pct.1 > 0.1 &
                               all_markers$avg_log2FC > 0.25 &
                               all_markers$p_val_adj < 0.05,][1:50,]
  
  result <- runCASSIA(
    marker_list = top_markers$gene,
    provider = "http://localhost:11434/v1",
    model = "gpt-oss:20b",
    tissue = "Brain",                 
    species = "Mouse" 
  )
  
  return(result)
})


