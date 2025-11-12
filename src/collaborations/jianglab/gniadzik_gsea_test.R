library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(ggplot2) # plotting
library(cowplot) # combining plots
library(data.table) # manipulate fgsea output

# for graphing and plotting
library(igraph)
library(ggraph)

# output directory
out_dir <- "~/Documents/projects/jianglab/weronika_gniadzik/gsea_results/"
dir.create(out_dir, showWarnings = F)

# Plotting parameters ----------------------------------------------------------
emap_pathways <- 20

cnet_pathways <- 10
cnet_genes <- 10


# Read in the files ------------------------------------------------------------
in_dir <- "~/Documents/projects/jianglab/weronika_gniadzik/"

# get all CSV files in the directory
csv_files <- list.files(in_dir, ".csv", full.names=T)

# read in and clean up data
data_list <- lapply(csv_files, function(file) {
  
  data <- read.csv(file, sep = ";", dec = ",")
  colnames(data) <- c("gene", "log2FC", "neg_logFDR")
  
  # value needed for GSEA ranked test
  data$signed_neg_logFDR <- data$neg_logFDR *
    sign(data$log2FC)
  
  # remove duplicates
  data <- data[order(data$neg_logFDR, decreasing=T),]
  data <- data[!duplicated(data$gene),]
  
  return(data)
})
names(data_list) <- gsub(".csv", "", basename(csv_files))

# Get pathways -----------------------------------------------------------------

reactome_gene_sets <- msigdbr(species = "human", collection = "C2", subcollection = "CP:REACTOME")
gobp_gene_sets <- msigdbr(species = "human", collection = "C5", subcollection = "GO:BP")

# clean up names
reactome_gene_sets$gs_name <- tolower(gsub("\\_", " ", gsub("REACTOME_", "", reactome_gene_sets$gs_name)))
gobp_gene_sets$gs_name <- tolower(gsub("\\_", " ", gsub("GOBP_", "", gobp_gene_sets$gs_name)))

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(reactome=list_convert(reactome_gene_sets),
                        gobp=list_convert(gobp_gene_sets))

# Perform GSEA Test ------------------------------------------------------------
gsea_results <- lapply(data_list, function(data) {
  
  # rank genes based on the signed negative logp
  ranked_genes <- data %>%
    filter(!is.na(signed_neg_logFDR)) %>%
    arrange(desc(signed_neg_logFDR)) %>%
    pull(signed_neg_logFDR, name = gene)
  
  # determine the GSEA results for each pathway set
  total_gsea_results <- lapply(names(total_gene_sets), function(gs_name) {
    
    print(gs_name)
    
    gene_sets <- total_gene_sets[[gs_name]]
    # run GSEA
    fgsea_results <- fgsea(
      pathways = gene_sets,
      stats = ranked_genes,
      minSize = 15, # Minimal size of a gene set to test
      maxSize = 500 # Maximal size of a gene set to test
    )
    
    # use the "leading edge" genes to determine an average fold change
    fgsea_results$avg_logfc <- lapply(fgsea_results$leadingEdge, function(gene_list) {
      
      mean(data[data$gene %in% gene_list,]$log2FC)
      
    })
    
    # reorder columns
    fgsea_results <- setcolorder(fgsea_results,
                                 c(colnames(fgsea_results)[colnames(fgsea_results) != "leadingEdge"],
                                     "leadingEdge"))
    
    return(fgsea_results[order(fgsea_results$pval),])
  })
  names(total_gsea_results) <- names(total_gene_sets)
  
  # return GSEA results as well as ranked list
  return(list(ranked_genes=ranked_genes,
              gsea_results=total_gsea_results))
})

# Write results to excel files
for (data_name in names(gsea_results)) {
  
  gsea_data <- gsea_results[[data_name]]
  
  write.xlsx(gsea_data$gsea_results,
             file=paste0(out_dir, data_name, ".gsea_results.xlsx"),
             colWidths="auto")
  
}

# EMAP Plots -------------------------------------------------------------------

# function to create the plots
create_fgsea_emapplot <- function(fgsea_res, 
                                  pathways_list, 
                                  top_n = 20,
                                  similarity_cutoff = 0.2,
                                  use_leading_edge = TRUE) {
  # Select top pathways
  top_pathways <- head(fgsea_res, top_n)
  pathway_names <- top_pathways$pathway
  
  # Calculate Jaccard similarity between pathways
  n_pathways <- length(pathway_names)
  similarity_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
  rownames(similarity_matrix) <- pathway_names
  colnames(similarity_matrix) <- pathway_names
  
  for (i in 1:n_pathways) {
    for (j in i:n_pathways) {
      # Use leading edge genes if specified, otherwise use full pathway
      if (use_leading_edge) {
        genes_i <- unlist(top_pathways$leadingEdge[i])
        genes_j <- unlist(top_pathways$leadingEdge[j])
      } else {
        genes_i <- pathways_list[[pathway_names[i]]]
        genes_j <- pathways_list[[pathway_names[j]]]
      }
      
      intersection <- length(intersect(genes_i, genes_j))
      union <- length(union(genes_i, genes_j))
      
      jaccard <- intersection / union
      similarity_matrix[i, j] <- jaccard
      similarity_matrix[j, i] <- jaccard
    }
  }
  
  # Create edge list for pathways with similarity > cutoff
  edges <- data.frame()
  for (i in 1:(n_pathways-1)) {
    for (j in (i+1):n_pathways) {
      if (similarity_matrix[i, j] > similarity_cutoff) {
        edges <- rbind(edges, data.frame(
          from = pathway_names[i],
          to = pathway_names[j],
          similarity = similarity_matrix[i, j]
        ))
      }
    }
  }
  
  # Create graph
  if (nrow(edges) == 0) {
    warning("No pathway connections above similarity cutoff. Try lowering similarity_cutoff.")
    # Create disconnected graph
    g <- graph_from_data_frame(
      data.frame(from = character(0), to = character(0)),
      directed = FALSE,
      vertices = pathway_names
    )
  } else {
    g <- graph_from_data_frame(edges, directed = FALSE, vertices = pathway_names)
    E(g)$weight <- edges$similarity
  }
  
  # Add node attributes
  V(g)$padj <- top_pathways$padj[match(V(g)$name, top_pathways$pathway)]
  V(g)$NES <- top_pathways$NES[match(V(g)$name, top_pathways$pathway)]
  V(g)$size <- top_pathways$size[match(V(g)$name, top_pathways$pathway)]
  V(g)$avg_logfc <- as.numeric(top_pathways$avg_logfc[match(V(g)$name, top_pathways$pathway)])
  
  # Plot
  set.seed(123)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = similarity, alpha = similarity), 
                   color = "grey60") +
    geom_node_point(aes(size = size, color = avg_logfc)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, 
                   max.overlaps = 20) +
    scale_edge_width_continuous(range = c(0.5, 2), name = "Similarity") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    scale_size_continuous(range = c(3, 10), name = "Gene Set Size") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "Avg Log2FC") +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    ggtitle("Enrichment Map of Top Pathways")
  
  return(list(plot = p, similarity_matrix = similarity_matrix))
}


# create directory for plots
emap_plot_dir <- paste0(out_dir, "emap_plots/")
dir.create(emap_plot_dir, showWarnings = F)

# loop through results and make plots
for (data_type in names(gsea_results)) {
  
  # get the GSEA results for the comparison
  data_gsea_results <- gsea_results[[data_type]]
  
  # create a plot for each geneset type
  for (gs_name in names(data_gsea_results$gsea_results)) {
    
    # get pathways and GSEA results
    pathways <- total_gene_sets[[gs_name]]
    fgsea_results <- data_gsea_results$gsea_results[[gs_name]]
    
    # make the plot
    emap <- create_fgsea_emapplot(fgsea_results,
                                  pathways,
                                  top_n = emap_pathways)
    
    # output and save the plot
    print(emap$plot)
    ggsave(paste0(emap_plot_dir, data_type, ".", gs_name, ".emap_pathways.png"),
                  width=8, height=6, bg = "white")
    
  }
  
}

# CNET plots ------------------------------------------------------------------

# Function to extract leading edge genes for significant pathways
prepare_cnetplot_data <- function(fgsea_res, 
                                  gene_data, 
                                  pathways_list, 
                                  padj_cutoff = 0.05) {
  # Filter significant results
  sig <- fgsea_res[padj < padj_cutoff]
  
  if (nrow(sig) == 0) {
    stop("No significant pathways found. Try increasing padj_cutoff.")
  }
  
  # Create gene-pathway mapping
  gene_pathway_list <- list()
  
  for (i in 1:nrow(sig)) {
    pathway_name <- sig$pathway[i]
    leading_genes <- unlist(sig$leadingEdge[i])
    gene_pathway_list[[pathway_name]] <- leading_genes
  }
  
  return(list(
    pathways = sig,
    gene_pathway_map = gene_pathway_list
  ))
}

# create the CNET plot
create_fgsea_cnetplot <- function(fgsea_res, gene_pathway_map, 
                                  gene_data, top_n = 5, 
                                  genes_per_pathway = 10) {
  # Select top pathways
  top_pathways <- head(fgsea_res, top_n)
  
  # Create edge list
  edges <- data.frame()
  
  for (pathway in top_pathways$pathway) {
    genes <- gene_pathway_map[[pathway]]
    # Limit genes per pathway for visualization
    genes <- head(genes, genes_per_pathway)
    
    for (gene in genes) {
      edges <- rbind(edges, data.frame(
        from = pathway,
        to = gene
      ))
    }
  }
  
  # Create graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Add node attributes
  V(g)$type <- ifelse(V(g)$name %in% top_pathways$pathway, "pathway", "gene")
  
  # Add fold change for genes
  V(g)$log2FC <- NA
  gene_idx <- which(V(g)$type == "gene")
  for (i in gene_idx) {
    gene_name <- V(g)$name[i]
    if (gene_name %in% gene_data$gene) {
      V(g)$log2FC[i] <- as.numeric(gene_data$log2FC[gene_data$gene == gene_name])
    }
  }
  
  # Plot
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.3, color = "grey70") +
    geom_node_point(aes(color = log2FC, size = ifelse(type == "pathway", 5, 3),
                        shape = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, na.value = "grey80",
                          name = "Fold Change") +
    scale_size_identity() +
    scale_shape_manual(values = c("pathway" = 15, "gene" = 19),
                       name = "Node Type") +
    theme_void() +
    theme(legend.position = "right") +
    ggtitle("Gene-Pathway Network")
  
  return(p)
}

# create directory for CNET plots
cnet_plot_dir <- paste0(out_dir, "cnet_plots/")
dir.create(cnet_plot_dir, showWarnings = F)

# loop through results and make plots
for (data_type in names(gsea_results)) {
  
  # get the GSEA results and the gene data for the comparison
  data_gsea_results <- gsea_results[[data_type]]
  gene_data <- data_list[[data_type]]
  
  # create a plot for each geneset type
  for (gs_name in names(data_gsea_results$gsea_results)) {
    # get pathways and GSEA results
    pathways <- total_gene_sets[[gs_name]]
    fgsea_results <- data_gsea_results$gsea_results[[gs_name]]
    
    # prepare the data
    cnet_data <- prepare_cnetplot_data(fgsea_results,
                                       gene_data,
                                       pathways)
    
    # make the plot
    create_fgsea_cnetplot(cnet_data$pathways,
                          cnet_data$gene_pathway_map,
                          gene_data,
                          top_n = cnet_pathways,
                          genes_per_pathway = cnet_genes)
    ggsave(paste0(cnet_plot_dir, data_type, ".", gs_name, ".cnet_plot.png"),
           width=8, height=6, bg = "white")
    
  }
}



