% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_clustermap.R
\name{return_score}
\alias{return_score}
\title{Retrieve importance score}
\usage{
return_score(df, motif_name, condition)
}
\arguments{
\item{df}{the data.frame coming from \code{models2dataframe}, containing the top
k motif in each condition and their importance score}

\item{motif_name}{the name of a single motif}

\item{condition}{the name of a single condition}
}
\value{
the importance score for that particular motif, in that particular
condition. If the score cannot be retrieved (i.e. the motif is not among the top k
motifs in that condition) 0 is assigned
}
\description{
This function takes as input a data.frame coming
from \code{models2dataframe} call, a single motif name and a single condition,
retrieving the importance score for that motif in that condition.
If the score cannot be retrieved (i.e. the motif is not among the top k
motifs in that condition) 0 is assigned
}
