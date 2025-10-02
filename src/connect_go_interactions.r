library(ontologyIndex)
library(openxlsx)
library(org.Hs.eg.db)
library(GO.db)
library(httr)
library(jsonlite)
library(biomaRt)

root_go_ids <- c("GO:0042552", # myelination
                 "GO:0031641", # regulation of myelination
                 "GO:0032288", # myelin assembly
                 "GO:0043217") # myelin maintenance

# get GO ontology
go <- get_ontology("https://purl.obolibrary.org/obo/go.obo",
                   extract_tags = "everything")

# get children of root GO terms
go_children <- get_descendants(go, root_go_ids)

# create Xref with IDs and term names
go_xref <- data.frame(go_id=go_children,
                      go_term=go$name[go_children])

# Connect to Ensembl
results <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_xref$go_id,
  columns = c("SYMBOL", "ENTREZID", "GENENAME", "EVIDENCE", "ONTOLOGY"),
  keytype = "GO"
)

go_results <- merge(go_xref,
                    results,
                    by.x="go_id",
                    by.y="GO")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

swissprot_ids <- getBM(
  attributes = c('hgnc_symbol', 'uniprotswissprot'),
  filters = 'hgnc_symbol',
  values = results$SYMBOL,
  mart = mart
)
swissprot_ids <- swissprot_ids[swissprot_ids$uniprotswissprot != "",]

# merge em in
go_results <- merge(go_results,
                    swissprot_ids,
                    by.x="SYMBOL",
                    by.y="hgnc_symbol")

# Common Evidence Codes
# 
# EXP: Inferred from Experiment
# IDA: Inferred from Direct Assay
# IPI: Inferred from Physical Interaction
# IMP: Inferred from Mutant Phenotype
# IGI: Inferred from Genetic Interaction
# IEP: Inferred from Expression Pattern
# TAS: Traceable Author Statement
# IC: Inferred by Curator
# NAS: Non-traceable Author Statement
# IEA: Inferred from Electronic Annotation (computationally assigned)

go_results$confidence <- ifelse(go_results$EVIDENCE %in% 
                                  c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"),
                                "high",
                                "low")

intact_results <- read.delim("data/stat4_intact_human.tsv")
# convert IDs
intact_results$id1 <- sapply(intact_results$X..ID.s..interactor.A, function(x) {
  unlist(strsplit(x, ":",))[2]
})
intact_results$id2 <- sapply(intact_results$ID.s..interactor.B, function(x) {
  unlist(strsplit(x, ":",))[2]
})

total_ids <- unique(c(intact_results$id1,
                      intact_results$id2))

# connected to Stat4
go_results_connect <- go_results[go_results$uniprotswissprot %in%
                                   total_ids,]

intact_results[intact_results$id1 == "P31749",]


