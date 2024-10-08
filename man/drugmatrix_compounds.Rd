% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxicogenomics-data.R
\name{drugmatrix_compounds}
\alias{drugmatrix_compounds}
\title{List of Compounds Available in DrugMatrix Database}
\usage{
drugmatrix_compounds(tissue = c("Liver", "Kidney", "Heart", "Hepatocytes"))
}
\arguments{
\item{tissue}{Optional argument specifying the tissue type. Default is missing. See the section
\strong{Available Data for Different Tissues} for details.

One of:
\itemize{
\item \code{"Liver"}, for data from liver tissue.
\item \code{"Kidney"}, for data from kidney tissue.
\item \code{"Heart"}, for data from heart tissue.
\item \code{"Hepatocytes"}, for data from rat primary hepatocytes.
}}
}
\value{
A data frame of boolean values if the parameter \code{tissue} is missing. The first column
contains available compound names. The subsequent columns represent different types of
available data. In the data frame:
\itemize{
\item \code{TRUE} indicates that data is available for the corresponding compound.
\item \code{FALSE} signifies that it is not.
For a given \code{tissue}, the output will be a vector containing the available compound names.
}
}
\description{
The \code{drugmatrix_compounds()} function retrieves a list of available compounds from
the \href{https://ntp.niehs.nih.gov/data/drugmatrix}{DrugMatrix} database. These compound names
can be used to download and process Affymetrix CEL images and metadata from
\href{https://www.ncbi.nlm.nih.gov/geo/}{GEO}.
}
\section{Available Data for Different Tissues}{

The \href{https://ntp.niehs.nih.gov/data/drugmatrix}{DrugMatrix} database includes tests on
over 600 distinct compounds, conducted in both \verb{in vivo} and \verb{in vitro}
experimental setups. Male Sprague Dawley rats were used in both \verb{in vivo} and \verb{in vitro}
experiments, with \verb{in vitro} studies specifically employing rat primary hepatocytes.
The available data for different tissues are:
\itemize{
\item Liver
\item Kidney
\item Heart
\item Hepatocytes (Rat Primary Hepatocytes)
}
}

\examples{
# Data available for compounds in `Liver`
drugmatrix_compounds(tissue = "Liver")
# Data available for compounds in `Kidney`
drugmatrix_compounds(tissue = "Kidney")
# Data available for compounds in `Heart`
drugmatrix_compounds(tissue = "Heart")
# Data available for compounds in `Hepatocytes`
drugmatrix_compounds(tissue = "Hepatocytes")
}
\seealso{
Other DrugMatrix data download helpers: 
\code{\link{drugmatrix_link}()},
\code{\link{drugmatrix_rawdata}()}
}
\concept{DrugMatrix data download helpers}
