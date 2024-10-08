% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppi-network.R
\name{get_hubdata}
\alias{get_hubdata}
\title{Extract hub network data}
\usage{
get_hubdata(net_data, condition = "degree >= 10", error_call = caller_env())
}
\arguments{
\item{net_data}{A network data of class \code{tgxtool}.}

\item{condition}{A character string specify the condition of hub network.
(default is \code{degree >= 20}). There are four metrics \code{degree}, \code{betweenness}, \code{closenes},
and \code{eigenes} can be used in the condition with the five comparison operators \code{<}, \code{>}, \code{<=},
\code{>=}, and \code{==}. At the end of the condition you must provide appropriate numeric value of metric
used in the condition.}
}
\value{
A network data
}
\description{
The function \code{get_hubdata()} is used to extract hub network data from full network data.
}
\examples{
\dontrun{
net_data <- get_netdata(
  gene_data,
  gene_col = "gene_name",
  p.value_col = "group",
  organism = "rat",
  score_threshold = 400,
  cluster_method = "edge.betweenness"
)
get_hubdata(net_data = net_dt, condition = "degree >= 20")
}
}
