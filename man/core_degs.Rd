% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core-genes.R
\name{core_degs}
\alias{core_degs}
\title{Identify Core Differentially Expressed Genes (DEGs)}
\usage{
core_degs(
  gene_df = res4,
  path_n = NULL,
  enrich_cutoff = 0.5,
  category = c("kegg", "go", "reactome", "wikipathways"),
  organism = c("rat", "human"),
  score_threshold = 200,
  version = "12"
)
}
\arguments{
\item{gene_df}{A data frame containing the gene information with a column named \code{gene_symbol}.}

\item{path_n}{Number of pathways to consider for enrichment analysis. Default is \code{NULL}.}

\item{enrich_cutoff}{Cutoff value for enrichment analysis. Default is \code{0.5}.}

\item{category}{Vector of categories for pathway enrichment analysis. Options are
\code{kegg}, \code{go}, \code{reactome}, and \code{wikipathways}. Default is \code{kegg} if not supplied.}

\item{organism}{Organism of interest. Options are \code{"rat"} and \code{"human"}. Default is \code{c("rat", "human")}.}

\item{score_threshold}{Score threshold for filtering genes. Default is \code{200}.}

\item{version}{Version of the database to use for pathway enrichment analysis. Default is \code{"12"}.}
}
\value{
A vector of pivot genes that are core DEGs across different clustering methods and enrichment analyses.
}
\description{
This function identifies core differentially expressed genes (DEGs) using
multiple clustering methods in PPI network and pathway enrichment analysis.
}
\details{
This function performs the following steps:
\itemize{
\item Validates the input categories and organisms.
\item Retrieves network data for the specified gene expression data using three clustering methods: "fastgreedy", "walktrap", and "edge.betweenness".
\item Performs pathway enrichment analysis for the specified categories and organism.
\item Identifies clusters in the network data and performs enrichment analysis on each cluster.
\item Calculates the proportion of enriched pathways shared with the overall enrichment analysis.
\item Filters clusters based on the enrichment cutoff and extracts unique genes from the filtered clusters.
\item Identifies and returns pivot genes that are common across the filtered clusters from different clustering methods.
}
}
\examples{
\dontrun{
  core_genes <- core_degs(gene_df = res4, path_n = 10, enrich_cutoff = 0.5, category = "kegg", organism = "human", score_threshold = 200, version = "12")
}

}
