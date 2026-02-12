library(Biostrings)
library(dplyr)

# read in fasta file
fasta_file <- "~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/Mus_musculus.GRCm39.cdna.all.fa.gz"

reads <- readDNAStringSet(fasta_file)

metadata <- lapply(names(reads), function(data) {
  
  split <- unlist(strsplit(data, " "))
  
  if (length(split) == 6) {
    
    return(data.frame(transcript_id=split[1],
                      type=split[2],
                      location=split[3],
                      gene_id=split[4],
                      gene_biotype=split[5],
                      transcript_biotype=split[6]))
    
  }
  
  if (!grepl("gene", split[7])) {
    
    return(data.frame(transcript_id=split[1],
                      type=split[2],
                      location=split[3],
                      gene_id=split[4],
                      gene_biotype=split[5],
                      transcript_biotype=split[6],
                      description=paste(split[7:length(split)], collapse=" ")))
    
  }
  
  data.frame(transcript_id=split[1],
             type=split[2],
             location=split[3],
             gene_id=split[4],
             gene_biotype=split[5],
             transcript_biotype=split[6],
             gene_symbol=split[7],
             description=paste(split[8:length(split)], collapse=" "))
  
})
metadata <- bind_rows(metadata)

# clean up the contents

metadata$gene_id <- gsub("gene:", "", metadata$gene_id)
metadata$gene_biotype <- gsub("gene_biotype:", "", metadata$gene_biotype)
metadata$transcript_biotype <- gsub("transcript_biotype:", "", metadata$transcript_biotype)
metadata$gene_symbol <- gsub("gene_symbol:", "", metadata$gene_symbol)
metadata$description <- gsub("description:", "", metadata$description)

# save results
write.csv(metadata,
          "~/Documents/projects/gaultierlab/sam_wachamo/bulkRNASeq/Mus_musculus.GRCm39.transcript2gene.csv",
          row.names = F)