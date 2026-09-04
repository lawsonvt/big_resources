library(openxlsx)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(snakecase)

root_dir <- "~/Documents/projects/ferrislab/cheyenne_harvey_stolz/"

out_dir <- paste0(root_dir, "results/data_explore/")
dir.create(out_dir, showWarnings = F, recursive = T)

# read in data
metadata <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Random number.xlsx"))
colnames(metadata)[1] <- "baf"
rownames(metadata) <- metadata$baf

n_matrix <- read.xlsx(paste0(root_dir, "BAF_CH_4681_Data.xlsx"), 
                      sheet = "p-value", startRow = 5, cols=1:51)

# add in sample type for ease

metadata$sample_type <- apply(metadata, 1, function(r) {
  
  if(grepl("Control 1", r[[4]])) {
    paste0(r[[3]], "c1")
  } else if (grepl("Control 2", r[[4]])) {
    paste0(r[[3]], "c2")
  } else if (grepl("Experimental 1", r[[4]])) {
    paste0(r[[3]], "e1")
  } else {
    paste0(r[[3]], "e2")
  }
  
})

# reformat n matrix

metabolite_xref <- n_matrix[,1:2]

metabolite_xref$metabolite_id <- 1:nrow(metabolite_xref)

dupes <- metabolite_xref[duplicated(metabolite_xref$Metabolite.name),]$Metabolite.name

# there are duplicates, so I can't simplify the matrix by using rownames
# WHY ARE THERE DUPLICATE METABOLITE NAMES?!?!


rownames(n_matrix) <- metabolite_xref$metabolite_id
rownames(metabolite_xref) <- metabolite_xref$metabolite_id

n_matrix$Metabolite.name <- NULL
n_matrix$Ontology <- NULL

# fix column names to match metadata baf
colnames(n_matrix) <- sapply(colnames(n_matrix), function(x) {
  as.numeric(unlist(strsplit(x, "_"))[1])
})

metadata[!metadata$baf %in% colnames(n_matrix),]

# remove from metadata missing sample
metadata <- metadata[metadata$baf %in% colnames(n_matrix),]

# do fold change analysis

#View(unique(metadata[,c("Group-detail","Internal.Treatment","sample_type")]))

contrasts <- c("E4e1-E4c1",
               "E4e2-E4c2",
               "E3e1-E3c1",
               "E2e1-E2c1",
               "KOe1-KOc1")

results_list <- lapply(contrasts, function(contrast) {
  
  print(contrast)
  
  group1 <- unlist(strsplit(contrast, "-"))[1]
  group2 <- unlist(strsplit(contrast, "-"))[2]
  
  group1_n <- n_matrix[,rownames(metadata[metadata$sample_type == group1,])]
  group2_n <- n_matrix[,rownames(metadata[metadata$sample_type == group2,])]
  
  group1_avg <- rowMeans(group1_n)
  group2_avg <- rowMeans(group2_n)
  
  log2fc <- log2(group1_avg / group2_avg)
  
  pvals <- sapply(rownames(n_matrix), function(i) {
    
    if (length(unique(c(group1_n[i,],group2_n[i,]))) == 1) {
      return(NA)
    }
    
    t.test(group1_n[i,],
           group2_n[i,],
           alternative = "two.sided",
           var.equal = T)$p.value
    
  })
  
  adj_pvals <- p.adjust(pvals, method="fdr")
  
  # output the results
  out <- data.frame(metabolite_xref[rownames(n_matrix),],
             contrast = contrast,
             group1_avg=group1_avg,
             group2_avg=group2_avg,
             log2fc=log2fc,
             pvalue=pvals,
             adj_pvalue=adj_pvals)
  out[order(out$pvalue),]
})
names(results_list) <- contrasts

write.xlsx(results_list, file=paste0(out_dir, "comparison_results.xlsx"),
           colWidths="auto")

write.xlsx(metadata, file=paste0(out_dir, "metadata.xlsx"))

# volcano plots

# one for each contrast

for (results in results_list) {
  
  contrast <- unique(results$contrast)
  
  results$log_p <- -log10(results$pvalue)
  
  sig_results <- results[results$pvalue < 0.05 &
                         abs(results$log2fc) > 0.5,]
  
  ggplot(results,
         aes(x=log2fc,
             y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = -log10(0.05),
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=sig_results,
               color="red",
               alpha=0.4) +
    geom_text_repel(data=sig_results,
                    aes(label=Metabolite.name),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", 
         title=contrast)
  ggsave(paste0(out_dir, contrast, ".volcano.png"), width=8, height=6)
  
}

# volcano for each ontological group

for (results in results_list) {
  
  contrast <- unique(results$contrast)
  
  results$log_p <- -log10(results$pvalue)
  sig_results <- results[results$pvalue < 0.05 &
                           abs(results$log2fc) > 0.5,]
  
  # make a directory for the volcanoes
  volcano_dir <- paste0(out_dir, contrast, ".ontology_volcanoes/")
  dir.create(volcano_dir, showWarnings = F)
  
  ontologies <- unique(results$Ontology)
  
  for (ontology in ontologies) {
    
    subset <- results[results$Ontology == ontology,]
    subset_sig <- sig_results[sig_results$Ontology == ontology,]
    
    ggplot(subset,
           aes(x=log2fc,
               y=log_p)) +
      geom_point(alpha=0.4, color="black") +
      geom_hline(yintercept = -log10(0.05),
                 color="red", linetype=2) +
      geom_vline(xintercept = 0.5,
                 color="red", linetype=2) +
      geom_vline(xintercept = -0.5,
                 color="red", linetype=2) +
      geom_point(data=subset_sig,
                 color="red",
                 alpha=0.4) +
      geom_text_repel(data=subset_sig,
                      aes(label=Metabolite.name),
                      color="red", size=2.5,
                      max.overlaps = 50) +
      theme_bw() +
      labs(x="Log2 Fold Change", y="-log10(P-Value)", 
           title=contrast, subtitle=ontology)
    ggsave(paste0(volcano_dir, contrast, ".", ontology, ".volcano.png"), width=5, height=4)
    
  }
    
}



