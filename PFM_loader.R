#' Load Position Frequency Matrices (PFMs) from a File
#'
#' This function reads a file containing position frequency matrices (PFMs) and
#' extracts them into a list of data frames. Each PFM represents the nucleotide
#' frequencies ("A", "C", "G", "T") at each position for a given motif.
#' 
#' @param path A character string representing the path to the file containing 
#' PFMs.The file should follow the MEME motif format, with motifs identified by 
#' the keyword "MOTIF".
#'
#' @return A named list of data frames. Each data frame corresponds to the PFM
#' of a motif and contains the nucleotide frequencies ("A", "C", "G", "T"). 
#' The names of the list elements are the motif names extracted from the file.
#' @importFrom utils read.table
#' @export


PFM_loader <- function(path) {
  Tab <- read.table(path, skip = 9, fill = TRUE, stringsAsFactors = FALSE)
  
  motif_indices <- grep("MOTIF", Tab$V1)
  motifs <- list()
  
  for (i in seq_along(motif_indices)) {
    start <- motif_indices[i]
    end <- if (i < length(motif_indices)) motif_indices[i + 1] - 1 else 
      nrow(Tab)
    
    motif_name <- Tab[start, 2]
    motif_matrix <- Tab[(start + 2):(end - 1), 1:4]
    colnames(motif_matrix) <- c("A", "C", "G", "T")
    motifs[[motif_name]] <- as.data.frame(apply(motif_matrix, 2, as.numeric), 
                                          stringsAsFactors = FALSE)
  }
  return(motifs)
}
