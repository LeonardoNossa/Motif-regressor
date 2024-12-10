#Funzione per ottenere tabella geneID geneName contenente primi geni
#tsv = TUSet.tsv
#pro = 511145.protein.info.v12.0.txt

table_first_genes <- function(tsv,pro){
  TUSet <- read.table(tsv, 
                      sep = "\t", 
                      header = TRUE, 
                      quote = "",
                      comment.char = "#",
                      stringsAsFactors = FALSE)
  
  # Funzione per filtrare la colonna X5.tuGenes
  tuGenes_filter <- function(x) {
    if ("X5.tuGenes" %in% colnames(x)) {
      x$X5.tuGenes <- sapply(strsplit(as.character(x$X5.tuGenes), ";"), `[`, 1)
    } else {
      stop("La colonna 'X5.tuGenes' non esiste")
    }
    return(x)
  }
  
  
  TUSet <- tuGenes_filter(TUSet)
  tuGenes_filtered <- TUSet$X5.tuGenes
  
  protein_info <- read.table(pro, 
                             sep = "\t", 
                             header = FALSE, 
                             quote = "",
                             comment.char = "#",
                             stringsAsFactors = FALSE)
  
  
  match_genes <- function(tuGenes_vector, protein_info_data) {
   
    matched_genes <- protein_info_data[protein_info_data[, 2] %in% tuGenes_vector, ]
    
    #string_protein_id (colonna 1) e preferred_name (colonna 2)
    result <- matched_genes[, c(1, 2)]
    return(result)
  }
  
  result_table <- match_genes(tuGenes_filtered, protein_info)
  
  #rinominare prima colonna
  rename_first_column <- function(dataframe) {
    dataframe$V1 <- lapply(strsplit(dataframe$V1, "\\."), function(x)x[2])
    return(dataframe)
  }
  
  result_table<-rename_first_column(result_table)
  
  colnames(result_table) [c(1,2)] <- c('Gene_ID', 'Gene_Name')
  
  return(result_table)
}
  




























