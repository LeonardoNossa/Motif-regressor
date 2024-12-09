# Importare avg_E_coli
data <- read.table("avg_E_coli_v4_Build_6_exps466probes4297.tab", 
                   sep = "\t", 
                   header = TRUE, 
                   quote = "",
                   comment.char = "#",
                   stringsAsFactors = FALSE)

# Modifica della prima colonna
rename_first_column <- function(vettore) {
  data[[1]] <- lapply(strsplit(data[[1]], "_"), function(x)x[2])
  return(data)
  }


data<-rename_first_column(data)






