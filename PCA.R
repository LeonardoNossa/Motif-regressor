install.packages("FactoMineR")
library(FactoMineR)

#' principal_comps
#'
#' @param matrix 
#' @param explained_var 
#' @param axis 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
principal_comps <- function(matrix,explained_var = 95, axis = 2, scale = TRUE){
  
  if(!(is.matrix(matrix))){
    stop("Please provide as input a matrix!")
  }
  
  if (axis == 1) {
    matrix <- t(matrix)
  }
  
  pca <- FactoMineR::PCA(matrix, scale.unit = scale, ncp = ncol(matrix), graph = FALSE)
  idx <- which(pca$eig[,3]>=explained_var)[1]
  final_coords <- pca$ind$coord[,seq_len(idx)]
  if (axis == 1){
    return(t(final_coords))
  }
  return(final_coords)
}

load("./Compendium_coli_bnum.RData")
compendium_red <- principal_comps(FINAL_TPM_FILTFILT, explained_var = 95, axis = 1, scale = TRUE)

# 99% --> 111 components
# 95% --> 51 components
# 90% --> 29 components

# we choose 95% (51 components)









