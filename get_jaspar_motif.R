#' get Jaspar motif
#'
#' This function queries the JASPAR database to retrieve position frequency 
#' matrices (PFMs) 
#' for transcription factor binding motifs of a given species, identified by 
#' its taxonomy ID.
#'
#' @param tax_id A numeric or character vector representing the taxonomy ID of 
#' the species 
#'               for which to retrieve motifs (e.g., 9606 for Homo sapiens).
#'
#' @return A named list of data frames, where each data frame represents the 
#' PFM of a motif.
#'         Each PFM contains the nucleotide frequencies ("A", "C", "G", "T") 
#'         at each position.
#'         The names of the list elements correspond to the transcription 
#'         factor names concatenated with their motif IDs 
#'         (e.g., "TF_name_MA0001").
#'
#' @import httr
#' @import jsonlite
#' @import glue
#' @export

get_jaspar_motif <- function(tax_id) {
  
  motif_id = c()
  url <- glue::glue("https://jaspar.elixir.no/api/v1/species/{tax_id}/")
  query <- list(identifiers = tax_id)
  
  response <- GET(url, query = query)
  content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
  
  while (TRUE){
    
    content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
    results <- content$results
    urls <- results$url
    url <- content$`next`
    if (length(url) != 0){
      response <- GET(url, query = query)
      motif_id <- c(motif_id, urls)
    } else {
      motif_id <- c(motif_id, urls)
      break
    }
  }
  motif_id <- vapply(X = motif_id, FUN = function(x, split) 
    {return(unlist(strsplit(x, split))[7])} , 
                     FUN.VALUE = character(1), split = "/")
  names(motif_id) <- c()
  
  matrices <- list()
  genes <- c()
  
  for (motif in motif_id) {
    url <- paste0("https://jaspar.elixir.no/api/v1/matrix/",motif,"/")
    response <- GET(url)
    if (status_code(response) == 200) {
      data <- content(response, "text", encoding = 'UTF-8')
      
      json_data <- fromJSON(data)
      
      genes <- c(genes, json_data$name)
      
      df <- data.frame("A"=json_data$pfm$A, "C"=json_data$pfm$C,
                       "G"=json_data$pfm$G, "T"=json_data$pfm$`T`)
      matrices <- append(matrices, list(df))
      
    } else {
      print(paste("Errore nella richiesta:", status_code(response)))
    }
    
  }
  
  names_df <- paste0(genes, "_", motif_id)
  names(matrices) <- names_df
  return(matrices)
  
}