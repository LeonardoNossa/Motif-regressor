#' Convert Position Frequency Matrices(PFMs) to Position Weight Matrices(PWMs)
#'
#' This function takes a list of Position Frequency Matrices (PFMs) and converts
#' them into Position Weight Matrices (PWMs) using the given background 
#' nucleotide frequencies.
#'
#' @param PFMs A named list of data frames, where each data frame represents a 
#' Position Frequency Matrix (PFM) for a motif. Each PFM should contain 
#' nucleotide frequencies for "A", "C", "G", and "T" at each position.
#' @param background A numeric vector of length 4 representing the background
#'  nucleotide probabilities for "A", "C", "G", and "T". 
#'  The default is `c(0.25, 0.25, 0.25, 0.25)`. Values should sum to 1.
#'
#' @return A named list of data frames, where each data frame corresponds to a
#'  Position Weight Matrix (PWM) for a motif. Each PWM contains log2-transformed
#'  weights relative to the background probabilities.
#' @export

PFM2PWM <- function(PFMs, background = c(0.25, 0.25, 0.25, 0.25)) {
  PWMs <- list()
  for (name in names(PFMs)) {
    p <- as.matrix(PFMs[[name]]) + 0.01
    p <- p / rowSums(p)
    p <- log2(p / background)
    PWMs[[name]] <- as.data.frame(p)
  }
  return(PWMs)
}