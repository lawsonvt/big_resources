library(ggplot2)
library(ggrepel)
library(snakecase)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
library(reshape2)
library(openxlsx)

root_dir <- "~/Documents/projects/sheybanilab/stef_maslova/"

out_dir <- paste0(root_dir, "results/gsva_method_comp.salmon/")
dir.create(out_dir, showWarnings = F)

# load in data sets
total_results <- readRDS(paste0(root_dir, "results/gsea_workflow.salmon/total.gsea_results.RDS"))
direct_results <- readRDS(paste0(root_dir, "results/gsea_workflow.salmon.direct_comparison/total.gsea_results.RDS"))

# combine them!

combined_results <- lapply(names(total_results), function(comparison) {
  
  pathways <- names(total_results[[comparison]]$gsea_results)
  
  merged_pathways <- lapply(pathways, function(pathway) {
    
    total_pathway <- total_results[[comparison]]$gsea_results[[pathway]]
    direct_pathway <- direct_results[[comparison]]$gsea_results[[pathway]]
    
    # add in metadata
    total_pathway$pathway_db <- pathway
    total_pathway$comparison <- comparison
    
    direct_pathway$pathway_db <- pathway
    direct_pathway$comparison <- comparison
    
    # merge em
    merge(total_pathway,
          direct_pathway,
          by=c("pathway","pathway_db","comparison"),
          suffixes=c(".total",".direct"))
    
  })
  names(merged_pathways) <- pathways
  
  return(merged_pathways)
})
names(combined_results) <- names(total_results)

# bind them
combined_results_df <- bind_rows(lapply(names(combined_results), function(comparison) {
  
  bind_rows(combined_results[[comparison]])
  
}))


# pull out the top pathways

max_fdr <- 0.05

top_hallmark <- combined_results_df[combined_results_df$pathway_db == "hallmark" &
                                      (combined_results_df$padj.total < max_fdr |
                                         combined_results_df$padj.direct < max_fdr),]

top_hallmark_pathways <- unique(top_hallmark$pathway)

hallmark_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_hallmark_pathways,],
                            pathway ~ comparison, value.var="NES.direct")
hallmark_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_hallmark_pathways,],
                             pathway ~ comparison, value.var="padj.direct")
hallmark_direct_fdr_print <- hallmark_direct_fdr

hallmark_direct_fdr_print[hallmark_direct_fdr < max_fdr] <- "*"
hallmark_direct_fdr_print[hallmark_direct_fdr > max_fdr] <- ""

hallmark_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_hallmark_pathways,],
                             pathway ~ comparison, value.var="NES.total")
hallmark_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_hallmark_pathways,],
                             pathway ~ comparison, value.var="padj.total")
hallmark_total_fdr_print <- hallmark_total_fdr

hallmark_total_fdr_print[hallmark_total_fdr < max_fdr] <- "*"
hallmark_total_fdr_print[hallmark_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

hallmark_direct_nes <- hallmark_direct_nes[,column_order]
hallmark_direct_fdr_print <- hallmark_direct_fdr_print[,column_order]

h_direct <- Heatmap(hallmark_direct_nes,
        col=colorRamp2(c(-2,0,2), c("blue","white","red")),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(hallmark_direct_fdr_print[i,j], x, y)
        }, name="NES",
        width= unit(5, "cm"),
        column_title = "Direct Comparison",
        cluster_columns = F)

hallmark_total_nes <- hallmark_total_nes[,column_order]
hallmark_total_fdr_print <- hallmark_total_fdr_print[,column_order]

h_total <- Heatmap(hallmark_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(hallmark_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F)

pdf(paste0(out_dir, "top_hallmark_pathways.total_direct.pdf"), width=12, height=9)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# try same for reactome

top_reactome <- combined_results_df[combined_results_df$pathway_db == "reactome" &
                                      (combined_results_df$padj.total < max_fdr |
                                         combined_results_df$padj.direct < max_fdr),]

# get reactome immune pathways

immune_pathways <- read.xlsx("resources/reactome_immune_pathways.xlsx")

immune_pathways$msigdb_name <- paste0("REACTOME_", toupper(to_snake_case(immune_pathways$Pathway.Name)))

sum(immune_pathways$msigdb_name %in% combined_results_df$pathway)

found <- immune_pathways[immune_pathways$msigdb_name %in% combined_results_df$pathway,]
not_found <- immune_pathways[!immune_pathways$msigdb_name %in% combined_results_df$pathway,]

msigdb_reactome_pathways <- unique(grep("REACTOME", combined_results_df$pathway, value=T))

# filter top to immune pathways
top_reactome <- top_reactome[top_reactome$pathway %in% immune_pathways$msigdb_name,]


top_reactome_pathways <- unique(top_reactome$pathway)

reactome_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_reactome_pathways,],
                             pathway ~ comparison, value.var="NES.direct")
reactome_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_reactome_pathways,],
                             pathway ~ comparison, value.var="padj.direct")
reactome_direct_fdr_print <- reactome_direct_fdr

reactome_direct_fdr_print[reactome_direct_fdr < max_fdr] <- "*"
reactome_direct_fdr_print[reactome_direct_fdr > max_fdr] <- ""

reactome_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_reactome_pathways,],
                            pathway ~ comparison, value.var="NES.total")
reactome_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_reactome_pathways,],
                            pathway ~ comparison, value.var="padj.total")
reactome_total_fdr_print <- reactome_total_fdr

reactome_total_fdr_print[reactome_total_fdr < max_fdr] <- "*"
reactome_total_fdr_print[reactome_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

reactome_direct_nes <- reactome_direct_nes[,column_order]
reactome_direct_fdr_print <- reactome_direct_fdr_print[,column_order]

h_direct <- Heatmap(reactome_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(reactome_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

reactome_total_nes <- reactome_total_nes[,column_order]
reactome_total_fdr_print <- reactome_total_fdr_print[,column_order]

h_total <- Heatmap(reactome_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(reactome_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_immune_reactome_pathways.total_direct.pdf"), width=12, height=6)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# GO BP
top_gobp <- combined_results_df[combined_results_df$pathway_db == "gobp" &
                                      (combined_results_df$padj.total < max_fdr |
                                         combined_results_df$padj.direct < max_fdr),]



go_immune_pathways <- read.xlsx("resources/go_immune_response_biological_processes.xlsx")

go_immune_pathways$msigdb_name <- paste0("GOBP_", toupper(to_snake_case(go_immune_pathways$GO.Biological.Process)))

sum(go_immune_pathways$msigdb_name %in% combined_results_df$pathway)


# filter top to immune pathways
top_gobp <- top_gobp[top_gobp$pathway %in% go_immune_pathways$msigdb_name,]


top_gobp_pathways <- unique(top_gobp$pathway)

gobp_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_gobp_pathways,],
                             pathway ~ comparison, value.var="NES.direct")
gobp_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_gobp_pathways,],
                             pathway ~ comparison, value.var="padj.direct")
gobp_direct_fdr_print <- gobp_direct_fdr

gobp_direct_fdr_print[gobp_direct_fdr < max_fdr] <- "*"
gobp_direct_fdr_print[gobp_direct_fdr > max_fdr] <- ""

gobp_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_gobp_pathways,],
                            pathway ~ comparison, value.var="NES.total")
gobp_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_gobp_pathways,],
                            pathway ~ comparison, value.var="padj.total")
gobp_total_fdr_print <- gobp_total_fdr

gobp_total_fdr_print[gobp_total_fdr < max_fdr] <- "*"
gobp_total_fdr_print[gobp_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

gobp_direct_nes <- gobp_direct_nes[,column_order]
gobp_direct_fdr_print <- gobp_direct_fdr_print[,column_order]

h_direct <- Heatmap(gobp_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(gobp_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

gobp_total_nes <- gobp_total_nes[,column_order]
gobp_total_fdr_print <- gobp_total_fdr_print[,column_order]

h_total <- Heatmap(gobp_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(gobp_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_immune_gobp_pathways.total_direct.pdf"), width=12, height=6)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# KEGG

top_kegg <- combined_results_df[combined_results_df$pathway_db == "kegg" &
                                      (combined_results_df$padj.total < max_fdr |
                                         combined_results_df$padj.direct < max_fdr),]


top_kegg_pathways <- unique(top_kegg$pathway)

kegg_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_kegg_pathways,],
                         pathway ~ comparison, value.var="NES.direct")
kegg_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_kegg_pathways,],
                         pathway ~ comparison, value.var="padj.direct")
kegg_direct_fdr_print <- kegg_direct_fdr

kegg_direct_fdr_print[kegg_direct_fdr < max_fdr] <- "*"
kegg_direct_fdr_print[kegg_direct_fdr > max_fdr] <- ""

kegg_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_kegg_pathways,],
                        pathway ~ comparison, value.var="NES.total")
kegg_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_kegg_pathways,],
                        pathway ~ comparison, value.var="padj.total")
kegg_total_fdr_print <- kegg_total_fdr

kegg_total_fdr_print[kegg_total_fdr < max_fdr] <- "*"
kegg_total_fdr_print[kegg_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

kegg_direct_nes <- kegg_direct_nes[,column_order]
kegg_direct_fdr_print <- kegg_direct_fdr_print[,column_order]

h_direct <- Heatmap(kegg_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(kegg_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

kegg_total_nes <- kegg_total_nes[,column_order]
kegg_total_fdr_print <- kegg_total_fdr_print[,column_order]

h_total <- Heatmap(kegg_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(kegg_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_kegg_pathways.total_direct.pdf"), width=12, height=8)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# top rankings

top_num <- 10

top_num_reactome_direct <- lapply(names(direct_results), function(comparison) {
  
  return(direct_results[[comparison]]$gsea_results$reactome$pathway[1:top_num])
  
})
top_num_reactome_direct <- unique(unlist(top_num_reactome_direct))

top_num_reactome_total <- lapply(names(total_results), function(comparison) {
  
  return(total_results[[comparison]]$gsea_results$reactome$pathway[1:top_num])
  
})
top_num_reactome_total <- unique(unlist(top_num_reactome_total))

# merge em up
top_num_reactome <- unique(c(top_num_reactome_direct,
                             top_num_reactome_total))



reactome_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_num_reactome,],
                             pathway ~ comparison, value.var="NES.direct")
reactome_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_num_reactome,],
                             pathway ~ comparison, value.var="padj.direct")
reactome_direct_fdr_print <- reactome_direct_fdr

reactome_direct_fdr_print[reactome_direct_fdr < max_fdr] <- "*"
reactome_direct_fdr_print[reactome_direct_fdr > max_fdr] <- ""

reactome_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_num_reactome,],
                            pathway ~ comparison, value.var="NES.total")
reactome_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_num_reactome,],
                            pathway ~ comparison, value.var="padj.total")
reactome_total_fdr_print <- reactome_total_fdr

reactome_total_fdr_print[reactome_total_fdr < max_fdr] <- "*"
reactome_total_fdr_print[reactome_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

reactome_direct_nes <- reactome_direct_nes[,column_order]
reactome_direct_fdr_print <- reactome_direct_fdr_print[,column_order]

h_direct <- Heatmap(reactome_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(reactome_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

reactome_total_nes <- reactome_total_nes[,column_order]
reactome_total_fdr_print <- reactome_total_fdr_print[,column_order]

h_total <- Heatmap(reactome_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(reactome_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_num_reactome_pathways.total_direct.pdf"), width=12, height=8)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# GO BP

top_num <- 10

top_num_gobp_direct <- lapply(names(direct_results), function(comparison) {
  
  return(direct_results[[comparison]]$gsea_results$gobp$pathway[1:top_num])
  
})
top_num_gobp_direct <- unique(unlist(top_num_gobp_direct))

top_num_gobp_total <- lapply(names(total_results), function(comparison) {
  
  return(total_results[[comparison]]$gsea_results$gobp$pathway[1:top_num])
  
})
top_num_gobp_total <- unique(unlist(top_num_gobp_total))

# merge em up
top_num_gobp <- unique(c(top_num_gobp_direct,
                             top_num_gobp_total))



gobp_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_num_gobp,],
                             pathway ~ comparison, value.var="NES.direct")
gobp_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_num_gobp,],
                             pathway ~ comparison, value.var="padj.direct")
gobp_direct_fdr_print <- gobp_direct_fdr

gobp_direct_fdr_print[gobp_direct_fdr < max_fdr] <- "*"
gobp_direct_fdr_print[gobp_direct_fdr > max_fdr] <- ""

gobp_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_num_gobp,],
                            pathway ~ comparison, value.var="NES.total")
gobp_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_num_gobp,],
                            pathway ~ comparison, value.var="padj.total")
gobp_total_fdr_print <- gobp_total_fdr

gobp_total_fdr_print[gobp_total_fdr < max_fdr] <- "*"
gobp_total_fdr_print[gobp_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

gobp_direct_nes <- gobp_direct_nes[,column_order]
gobp_direct_fdr_print <- gobp_direct_fdr_print[,column_order]

h_direct <- Heatmap(gobp_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(gobp_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

gobp_total_nes <- gobp_total_nes[,column_order]
gobp_total_fdr_print <- gobp_total_fdr_print[,column_order]

h_total <- Heatmap(gobp_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(gobp_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_num_gobp_pathways.total_direct.pdf"), width=12, height=10)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()

# wiki pathways?


top_wiki <- combined_results_df[combined_results_df$pathway_db == "wikipathways" &
                                  (combined_results_df$padj.total < max_fdr |
                                     combined_results_df$padj.direct < max_fdr),]


top_wiki_pathways <- unique(top_wiki$pathway)

wiki_direct_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_wiki_pathways,],
                         pathway ~ comparison, value.var="NES.direct")
wiki_direct_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_wiki_pathways,],
                         pathway ~ comparison, value.var="padj.direct")
wiki_direct_fdr_print <- wiki_direct_fdr

wiki_direct_fdr_print[wiki_direct_fdr < max_fdr] <- "*"
wiki_direct_fdr_print[wiki_direct_fdr > max_fdr] <- ""

wiki_total_nes <- acast(combined_results_df[combined_results_df$pathway %in% top_wiki_pathways,],
                        pathway ~ comparison, value.var="NES.total")
wiki_total_fdr <- acast(combined_results_df[combined_results_df$pathway %in% top_wiki_pathways,],
                        pathway ~ comparison, value.var="padj.total")
wiki_total_fdr_print <- wiki_total_fdr

wiki_total_fdr_print[wiki_total_fdr < max_fdr] <- "*"
wiki_total_fdr_print[wiki_total_fdr > max_fdr] <- ""

column_order <- c("CAR-Control",
                  "FUS-Control",
                  "CAR_FUS-Control")

wiki_direct_nes <- wiki_direct_nes[,column_order]
wiki_direct_fdr_print <- wiki_direct_fdr_print[,column_order]

h_direct <- Heatmap(wiki_direct_nes,
                    col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                    cell_fun = function(j, i, x, y, w, h, col) {
                      grid.text(wiki_direct_fdr_print[i,j], x, y)
                    }, name="NES",
                    width= unit(5, "cm"),
                    column_title = "Direct Comparison",
                    cluster_columns = F,row_names_gp = gpar(fontsize = 8))

wiki_total_nes <- wiki_total_nes[,column_order]
wiki_total_fdr_print <- wiki_total_fdr_print[,column_order]

h_total <- Heatmap(wiki_total_nes,
                   col=colorRamp2(c(-2,0,2), c("blue","white","red")),
                   cell_fun = function(j, i, x, y, w, h, col) {
                     grid.text(wiki_total_fdr_print[i,j], x, y)
                   }, name="NES",
                   width = unit(5, "cm"),
                   column_title = "Total Samples",
                   cluster_columns = F,row_names_gp = gpar(fontsize = 8))

pdf(paste0(out_dir, "top_wiki_pathways.total_direct.pdf"), width=12, height=9)
draw(h_direct + h_total, heatmap_legend_side = "left")
dev.off()


