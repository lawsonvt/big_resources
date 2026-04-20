library(openxlsx)
library(tidyr)
library(prolfqua)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(snakecase)

# prolqua code adapted from:
# https://fgcz.github.io/prolfqua/index.html
# https://fgcz.github.io/prolfqua/articles/Comparing2Groups.html

# folder for analysis (you will need to change for your computer)
root_folder <- "~/Documents/projects/jianglab/ellie_sherman/"

# folder to store outputs
out_dir <- paste0(root_folder, "proteomics_analysis/")
dir.create(out_dir, showWarnings = F)

# read in the data
data <- read.xlsx(paste0(root_folder, "20260216_132640_APOE2 organoid_Report_Lulu Jiang lab.xlsx"),
                  startRow = 2)

# pivot data from wide format (current) to long format
data_long <- pivot_longer(data, cols = matches("(Quantity|Count)"), 
                          values_transform = as.numeric,
                          names_to="sample_raw")
# take NA counts out
na_counts <- data_long[data_long$sample_raw == "Count",]
data_long <- data_long[data_long$sample_raw != "Count",]


# fix sample ID
data_long$sample_id <- sapply(data_long$sample_raw,
                              function(x) {
                                
                                split_x <- unlist(strsplit(x, "[_.]"))
                                paste0(split_x[7:8], collapse="_")
                                
                              })

# make groups
data_long$group <- "APOE3_Ctrl"
data_long[grepl("A3AD", data_long$sample_id),]$group <- "APOE3_AD"
data_long[grepl("AZDU", data_long$sample_id),]$group <- "APOE2_Ctrl"
data_long[grepl("AZAD", data_long$sample_id),]$group <- "APOE2_AD"

# rename for analysis
data_long$protein_Id <- data_long$PG.ProteinAccessions

# reorg and simplify data
protein_xref <- unique(data_long[,c("protein_Id", "PG.Genes")])

data_long <- data_long[,c("sample_id","group","protein_Id", "value")]

# start prolfqua
# set configuration object
config <- prolfqua::AnalysisConfiguration$new()
config$file_name = "sample_id"
config$work_intensity = "value"
config$hierarchy[["protein_Id"]]    <-  "protein_Id"
config$factors[["group"]] <- "group"

# Build LFQData object
analysis_data <- prolfqua::setup_analysis(data_long, config)
lfqdata <- prolfqua::LFQData$new(analysis_data, config)

# make a plotter object
lfqplotter <- lfqdata$get_Plotter()

# Heatmap of NAs to see how they cluster
nah <- lfqplotter$na_heatmap()

pdf(paste0(out_dir, "na_heatmap.pdf"), width=7, height=6)
nah
dev.off()

# normalize and z-scale the data
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$config$is_response_transformed

# plotter for the transformed data
pl <- transformed$get_Plotter()

# detailed correlation plots between all samples
pl$pairs_smooth()

# overall heatmap correlation plot
p <- pl$heatmap_cor()

pdf(paste0(out_dir, "correlation_heatmap.pdf"), width=7, height=6)
p
dev.off()
# NOTE: THIS SHOWS A POTENTIAL OUTLIER!

# formula for modeling
formula_Condition <-  strategy_lm("transformedIntensity ~ group")

# specify model definition
modelName  <- "Model"
# Analysis contrasts
contr_spec <- c("APOE3_ADvsAPOE3_Ctrl" = "groupAPOE3_AD - groupAPOE3_Ctrl",
                "APOE2_ADvsAPOE2_Ctrl" = "groupAPOE2_AD - groupAPOE2_Ctrl")

# build model
mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_id = transformed$config$hierarchy_keys() )

# do contrasts
contr <- prolfqua::Contrasts$new(mod, contr_spec)

# create volcano plot (also runs the contrast analysis)
v1 <- contr$get_Plotter()$volcano()

# moderated Wald Test (recommended for simple contrasts like this one)
# also, it increases the amount of significant genes
contr <- prolfqua::ContrastsModerated$new(contr)

# pull out the modeled data
contrdf <- contr$get_contrasts()

# create another volcano plot
plotter <- contr$get_Plotter()
v2 <- plotter$volcano()

# plot comparison of volcano plots
gridExtra::grid.arrange(v1$FDR,v2$FDR, ncol = 1)

# plot and save comparison of volcano plots
plot_grid(v1$FDR,v2$FDR, ncol=1)
ggsave(paste0(out_dir, "volcanoes.png"), width=7, height=6)

# do it again, but drop the outlier
data_long <- data_long[data_long$sample_id != "A3AD_1",]

# Build LFQData object (use same config as before)
analysis_data <- prolfqua::setup_analysis(data_long, config)
lfqdata <- prolfqua::LFQData$new(analysis_data, config)

# make a plotter
lfqplotter <- lfqdata$get_Plotter()

# NA heatmap
nah <- lfqplotter$na_heatmap()

pdf(paste0(out_dir, "na_heatmap.drop_A3AD_1.pdf"), width=7, height=6)
nah
dev.off()

# normalize and scale data
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$config$is_response_transformed

# plotter for transformed data
pl <- transformed$get_Plotter()

# detailed correlation plots between all samples
pl$pairs_smooth()

# Correlation heatmap
p <- pl$heatmap_cor()

pdf(paste0(out_dir, "correlation_heatmap.drop_A3AD_1.pdf"), width=7, height=6)
p
dev.off()

# formula for model
formula_Condition <-  strategy_lm("transformedIntensity ~ group")

# specify model definition
modelName  <- "Model"
# Contrasts for analysis
contr_spec <- c("APOE3_ADvsAPOE3_Ctrl" = "groupAPOE3_AD - groupAPOE3_Ctrl",
                "APOE2_ADvsAPOE2_Ctrl" = "groupAPOE2_AD - groupAPOE2_Ctrl")

# build model
mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_id = transformed$config$hierarchy_keys() )

# create volcano plot (also runs the contrast analysis)
v1 <- contr$get_Plotter()$volcano()

# moderated Wald Test (recommended for simple contrasts like this one)
# also, it increases the amount of significant genes
contr <- prolfqua::ContrastsModerated$new(contr)

# create another volcano plot
plotter <- contr$get_Plotter()
v2 <- plotter$volcano()

# plot comparison of volcano plots
gridExtra::grid.arrange(v1$FDR,v2$FDR, ncol = 1)

# plot and save comparison of volcano plots
plot_grid(v1$FDR,v2$FDR, ncol=1)
ggsave(paste0(out_dir, "volcanoes.drop_A3AD_1.png"), width=7, height=6)


# get the results data
res <- contr$get_contrasts()

# merge back in gene names
res <- merge(res, protein_xref, by="protein_Id")

# clean up data
res <- res[,c("protein_Id","PG.Genes","diff","statistic","p.value","FDR","std.error","avgAbd","df",
              "conf.low","conf.high","sigma","contrast","modelName")]
colnames(res)[3] <- "log2fc"

# split for output 
res <- res[order(res$p.value),] # order results by P value

res_list <- list("APOE2_ADvsAPOE2_Ctrl"=res[res$contrast == "APOE2_ADvsAPOE2_Ctrl",],
                 "APOE3_ADvsAPOE3_Ctrl"=res[res$contrast == "APOE3_ADvsAPOE3_Ctrl",])

# Write to excel file
write.xlsx(res_list, file=paste0(out_dir,"contrast_results.xlsx"))


# Make better, prettier volcano plots

top_genes <- 25 # how many top gene names to show in plot

volcano_plot_list <-  lapply(names(res_list), function(contrast) {
  
  # loop through contrasts making volcano plots  
  subset <- res_list[[contrast]]
  
  subset$log_p <- -log10(subset$FDR)
  
  # significant proteins
  subset_sig <- subset[subset$FDR < 0.1 &
                         abs(subset$log2fc) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- 1
  
  p1 <- ggplot(subset,
               aes(x=log2fc,
                   y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = logp_thresh,
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=PG.Genes),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(FDR)", title=contrast)
  p1
  # plot individual volcanoes
  ggsave(paste0(out_dir, to_snake_case(contrast), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})

# one big plot of both volcanoes
plot_grid(plotlist = volcano_plot_list, ncol=2)
ggsave(paste0(out_dir, "combined_volcano_plots.png"), width=11, height=6)




