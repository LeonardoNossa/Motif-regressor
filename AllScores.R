#' Compute All PWM Scores
#'
#' This function calculates the maximum scores for a list of Position Weight 
#' Matrices (PWMs) across a set of DNA sequences using parallel computation. 
#' For each PWM, it computes the maximum score for every sequence and returns 
#' a matrix of scores.
#'
#' @param PWMs A list of Position Weight Matrices (PWMs). Each PWM in the list
#'  should be a numeric matrix with 4 rows (corresponding to A, C, G, T).
#' @param dataSeq A named character vector containing DNA sequences to be 
#' scored. Each sequence must consist of the characters "A", "C", "G", and "T". 
#' The names of the vector will be used as row names in the output.
#'
#' @return A numeric matrix where:
#' \item{Rows}{Correspond to the DNA sequences (named from `dataSeq`).}
#' \item{Columns}{Correspond to each PWM in `PWMs`. Each entry is the maximum 
#' score for the given sequence and PWM.}
#'
#' @import future
#' @importFrom future.apply future_lapply
#' @export
AllScores <- function(PWMs, dataSeq){
  plan(multisession, workers = 4) 
  Scores <- future_lapply(PWMs, function(pwm) scorer(pwm, dataSeq))
  Total <- do.call(cbind, lapply(Scores, function(score) score[, 2]))
  plan(sequential)
  rownames(Total) = names(dataSeq)
  return(Total)
}
