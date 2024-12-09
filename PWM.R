library(seqinr)
library(future)
library(future.apply)
library(glmnet)
library(e1071)
library(Biostrings)

path = "./motifregressor/jaspar.meme.txt"
path = "./motifregressor/SwissRegulon_e_coli.meme"

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