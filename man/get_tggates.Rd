% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxicogenomics-data.R
\name{get_tggates}
\alias{get_tggates}
\title{Retrieve and Process TG-GATES Data}
\usage{
get_tggates(
  compounds,
  species = c("Rat", "Human"),
  data_type = c("in_vivo", "in_vitro"),
  tissue = c("Liver", "Kidney"),
  dose_type = c("Single", "Repeat"),
  fc = TRUE,
  method = "rma",
  multicore = FALSE,
  store = FALSE,
  output_dir = missing_arg(),
  ...
)
}
\arguments{
\item{compounds}{A character vector of compound names to retrieve from the TG-GATES database.}

\item{data_type}{A character string specifying the type of data to retrieve. See \code{tggates_rawdata} for available options.}

\item{fc}{A logical value indicating whether to apply fold-change calculations. Default is \code{TRUE}.}

\item{method}{A character string specifying the preprocessing method to use. Default is \code{"rma"}.}

\item{multicore}{A logical value indicating whether to use multiple cores for processing. Default is \code{FALSE}.}

\item{store}{A logical value indicating whether to store the processed data. Default is \code{FALSE}.}

\item{output_dir}{A character string specifying the directory to store output files. Default is \code{missing_arg()}.}

\item{...}{Additional arguments passed to the \code{process_tggates} function.}
}
\value{
A list containing processed data from the TG-GATES database.
}
\description{
The \code{get_tggates()} function retrieves raw data for specified compounds from the TG-GATES database and processes it into gene expression values and metadata.
}
\examples{
compounds <- c("2-nitrofluorene", "N-nitrosomorpholine")
\dontrun{
data <- get_tggates(compounds, species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
head(data$expression)
head(data$metadata)
}
}
