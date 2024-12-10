#' Compute PWM Scores
#'
#' This function calculates the maximum position-specific scoring matrix (PSSM)
#' scores for a list of DNA sequences using a given Position Weight Matrix(PWM).
#' It identifies the position with the highest score for each sequence and 
#' returns both the position and the corresponding score.
#'
#' @param PWM A numeric matrix representing the Position Weight Matrix (PWM). 
#' The matrix should have 4 rows (for A, C, G, T).
#' @param dataSeq A character vector containing DNA sequences to be scored. 
#' Each sequence must be composed of the characters "A", "C", "G", and "T".
#'
#' @return A numeric matrix with two columns and a number of rows equal to the 
#' length of `dataSeq`:
#' \item{Column 1}{The position in the sequence with the highest score.}
#' \item{Column 2}{The maximum score for the sequence.}
#'
#' @export

scorer <- function(PWM,dataSeq){
  
  n_seq = length(dataSeq)
  Result <-matrix(0,ncol = 2,nrow= n_seq)
  PWM = matrix(as.numeric(as.matrix(PWM)),ncol = 4)
  
  for (s in (1: n_seq)) {
    UPSTREAM = s2c(dataSeq[s])
    MATRIX<-matrix(0,nrow=4,ncol=length(UPSTREAM))
    mA<-which(UPSTREAM=="A")
    mC<-which(UPSTREAM=="C")
    mG<-which(UPSTREAM=="G")
    mT<-which(UPSTREAM=="T")
    
    MATRIX[1,mA]<-1
    MATRIX[2,mC]<-1
    MATRIX[3,mG]<-1
    MATRIX[4,mT]<-1
    
    n_col = dim(PWM)[2] 
    ScoreMatrix<-matrix(0,ncol = 1, nrow=(length(UPSTREAM)-n_col))
    
    for (i in 1:(dim(MATRIX)[2] - n_col)-1){
      seq = MATRIX[,i:(nrow(PWM))]
      score = sum(diag(PWM %*% seq))
      ScoreMatrix[i,1] <- score
    }
    MaxIndex = which(ScoreMatrix == max(ScoreMatrix))
    Result[s,1] = MaxIndex[1]
    Result[s,2] = max(ScoreMatrix)[1]
  }
  return(Result)
}