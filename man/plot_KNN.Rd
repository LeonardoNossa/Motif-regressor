% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_regression.R
\name{plot_KNN}
\alias{plot_KNN}
\title{Plot KNN Model Metrics and Feature Importance}
\usage{
plot_KNN(KNN_model)
}
\arguments{
\item{KNN_model}{A trained KNN model object from the caret package.}
}
\value{
A list containing data frames and ggplot objects
for RMSE, R2, MAE, and feature importance plots.
}
\description{
Generates plots and data summaries for K-Nearest Neighbors (KNN)
model metrics (RMSE, R2, MAE) and feature importance.
}
