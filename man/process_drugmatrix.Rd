% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxicogenomics-data.R
\name{process_drugmatrix}
\alias{process_drugmatrix}
\title{Process Gene Expression and Metadata from DrugMatrix Database}
\usage{
process_drugmatrix(
  rawdata_path,
  fc = TRUE,
  method = "rma",
  multicore = FALSE,
  store = FALSE,
  output_dir = missing_arg(),
  ...
)
}
\arguments{
\item{rawdata_path}{A string specifying the full path to the raw data directory of compounds, which should contain a folder named \code{celfile} and a \code{tsv} attribute file.}

\item{fc}{A boolean indicating the format of the output gene expression data.
\itemize{
\item \code{TRUE} (default): The output will be in fold-change (FC) format.
\item \code{FALSE}: The output will be raw gene expression data.
}}

\item{method}{The normalization method to use, either \code{"rma"} or \code{"mas5"}. Default is \code{"rma"}.}

\item{multicore}{Used for parallel computing.
\itemize{
\item \code{FALSE} (default): No parallelization will be used.
\item \code{TRUE}: Computation will employ parallel processing, leaving one core for system use.
\item Integer: Specifies the number of workers for parallel processing.
}}

\item{store}{Logical. If \code{TRUE} (the default), the resulting data will be saved as
two CSV files (gene expression and metadata) in the output directory.}

\item{output_dir}{Output directory for downloading data. If not specified,
the temporary path of the current R session will be used.}

\item{...}{Additional arguments passed to \code{process_celfile}.}
}
\value{
A list with the following components:
\itemize{
\item \code{expression}: A matrix of gene expression values with probes in rows and samples (CEL IDs) in columns.
\item \code{metadata}: A data frame containing metadata about the samples.
}
}
\description{
The \code{process_drugmatrix()} function processes gene expression data and metadata for samples
from CEL files and attribute files in the DrugMatrix database.
}
\examples{
\dontrun{
file_path <- drugmatrix_rawdata(compounds = "Ethanol", "Liver")
proc_dm <- process_drugmatrix(file_path)
head(proc_dm$expression)
head(proc_dm$metadata)
}
}
\seealso{
Other Open TG-GATEs data download helpers: 
\code{\link{process_celfile}()},
\code{\link{process_tggates}()},
\code{\link{tggates_rawdata}()}
}
\concept{Open TG-GATEs data download helpers}
