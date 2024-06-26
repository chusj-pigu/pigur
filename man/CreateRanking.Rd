% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PathwayAnalysis.r
\name{CreateRanking}
\alias{CreateRanking}
\title{Create a Ranking of Genes from DE results.  The ranking is a based on a multiplication of the log2FC and the -log10(pvalue)}
\usage{
CreateRanking(
  deResultDF,
  pvalColumn = "p.value",
  log2FCColumn = "log2FC",
  geneColumn = "gene",
  twoSided = TRUE,
  groupByColumn = NULL
)
}
\arguments{
\item{deResultDF}{Data Frame containing the Differential expression results}

\item{pvalColumn}{Column name containing the p.value}

\item{log2FCColumn}{Column name containing the log2FC}

\item{geneColumn}{Column name containing the gene name}

\item{twoSided}{Is the ranking 2 sided (default) or not.  If not 2-sided, the absoluted value of the 2-sided ranking is returned.}

\item{groupByColumn}{If the input deResultDF contains multiple sets of results, the group by column should be specified to identify each unique resultset in the dataframe.  A ranking will be defined for each resultset.}
}
\value{
A list of rankings, one per resultset defined by the groupByColumn.  If the groupByColumn column is NULL, a single ranking is in the returned list and is named "All"
}
\description{
Create a Ranking of Genes from DE results.  The ranking is a based on a multiplication of the log2FC and the -log10(pvalue)
}
\examples{
\dontrun{
rankings = CreateRanking(myDEResultsDF,pvalColumn='p_val',log2FCColumn="avg_Log2FC",geneColumn='gene',groupByColumn="cluster")
names(rankings)
rankings[[1]]
}
}
