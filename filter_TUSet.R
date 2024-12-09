#Apro il file TUSet.tsv butto via tutti i geni che non sono la prima unità trascrizionale. 
#quindi guardo la colonna tugenes, dove ci sono più sequenze sono operoni;di questi devo tenere solo i primi se ce ne sono più di uno.
#poi match col compendium cosi anche li quali mi interessano


TUSet <- read.table("TUSet.tsv", 
                   sep = "\t", 
                   header = TRUE, 
                   comment.char = "#")


# Funzione per modificare la colonna 'X5.tuGenes'
tuGenes_filter <- function(TUSet) {

  if ("X5.tuGenes" %in% colnames(TUSet)) {
    # Modifica la colonna 'X5.tuGenes'
    TUSet$X5.tuGenes <- sapply(strsplit(as.character(TUSet$X5.tuGenes), ";"), `[`, 1)
  } else {
    stop("La colonna 'X5.tuGenes' non esiste")
  }
  return(TUSet)
}

TUSet <- tuGenes_filter(TUSet)

