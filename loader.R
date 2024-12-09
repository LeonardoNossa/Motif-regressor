library(seqinr)
library(future)
library(future.apply)
library(glmnet)
library(e1071)
library(Biostrings)

path = "./motifregressor/jaspar.meme.txt"
path = "./motifregressor/SwissRegulon_e_coli.meme"

PFM_loader <- function(path){
  
  Tab <-  read.table(path , skip = 9, fill = T)
  motif_indices <- grep("MOTIF", Tab$V1)
  if (all(nchar(Tab[Tab$V1 == "MOTIF",]$V3) > 0)) {
    motif_names <- paste0(Tab[Tab$V1 == "MOTIF",]$V3,"_",Tab[Tab$V1 == "MOTIF",]$V2)
  } else {
    motif_names <- Tab[Tab$V1 == "MOTIF",]$V2
  }
  

  motifs <- list()
  for (i in seq_along(motif_indices)) {
    
    start <- motif_indices[i]
    end <- if (i < length(motif_indices)) motif_indices[i + 1] - 1 else nrow(Tab)
    
    motif_name <- motif_names[i]
    motif_df <- Tab[(start + 2) : (end - 1),1:4]
    colnames(motif_df) <- c("A", "C", "G", "T")
    rownames(motif_df) <- seq(nrow(motif_df))
    
    if (any(motif_df == "URL")) {
      idx <- which(motif_df == "URL")
      motif_df <- motif_df[-(idx:nrow(motif_df)),]
    }
    
    motif_df <- as.data.frame(apply(motif_df,2,as.numeric), stringsAsFactors = FALSE)
    motifs[[motif_name]] <- motif_df
  }
  return(motifs)
}
