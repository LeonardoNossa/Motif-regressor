% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regressor.R
\name{Regression_KNN}
\alias{Regression_KNN}
\title{K-Nearest Neighbors (KNN) Regression}
\usage{
Regression_KNN(Scores, TPMs, Condition, params)
}
\arguments{
\item{Scores}{A matrix or data frame of scores to be used as predictors in
the KNN regression. This is internally converted to a matrix.}

\item{TPMs}{A data frame or matrix of TPM (Transcripts Per Million) values.
Columns correspond to different conditions or samples.}

\item{Condition}{A character string specifying the column name in \code{TPMs}
corresponding to the condition to be analyzed.}

\item{params}{A named list or vector containing parameters for the regression
model:
\describe{
\item{\code{num_cores}}{(Optional) Number of CPU cores to use for
parallelization. Default is the total number of available cores.
If a value greater than the available cores is provided, it will be
adjusted accordingly.}
\item{\code{cv}}{(Optional) Number of cross-validation folds. Default is 3.
If a value greater than the number of observations is provided, it will
default to Leave-One-Out Cross-Validation (LOOCV).}
}}
}
\value{
An object of class \code{RFE} representing the fitted KNN regression
model.
}
\description{
This function performs KNN regression using the \code{caret} package.
The best model is evaluated using Recursive Feature Elimination (RFE).
Computation of the best model is parallelized using \code{doParallel}.
}
