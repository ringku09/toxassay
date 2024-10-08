% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxicogenomics-data.R
\name{tggates_compounds}
\alias{tggates_compounds}
\title{List of Compounds Available in the TG-GATES Database}
\usage{
tggates_compounds(
  species = c("Rat", "Human"),
  data_type = c("in_vivo", "in_vitro"),
  tissue = c("Liver", "Kidney"),
  dose_type = c("Single", "Repeat"),
  error_call = caller_env()
)
}
\arguments{
\item{species}{A character vector specifying the species, either \code{"Rat"} or \code{"Human"}.}

\item{data_type}{A character vector specifying the type of data, either \code{"in_vivo"} or \code{"in_vitro"}.}

\item{tissue}{A character vector specifying the organ, either \code{"Liver"} or \code{"Kidney"}.}

\item{dose_type}{A character vector specifying the dose type for \emph{in-vivo} data, either \code{"Single"} or \code{"Repeat"}.}
}
\value{
A vector of compound names available for the specified criteria.
}
\description{
The \code{tggates_compounds()} function retrieves a list of available compounds from
the \href{https://toxico.nibiohn.go.jp/english/}{Open TG-GATEs} database. These compounds
can be used to download and process Affymetrix CEL images of samples.
}
\section{Open TG-GATEs Database}{

The \href{https://toxico.nibiohn.go.jp/english/}{TG-GATEs} database includes tests on 170
distinct compounds, primarily medicinal drugs, across both \emph{in-vivo} and \emph{in-vitro}
experimental setups. For \emph{in-vivo} trials, male Sprague Dawley rats were used, and
the studies were organized into single-dose and repeated-dose protocols. Additionally,
the database includes two \emph{in-vitro} studies utilizing primary hepatocytes from male
Sprague Dawley rats and human donors.
}

\examples{
# Data available for compounds in Rat in-vivo liver Single data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
# Data available for compounds in Rat in-vivo liver Repeat data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Repeat")
# Data available for compounds in Rat in-vitro liver Single data
tggates_compounds(species = "Rat", data_type = "in_vitro", tissue = "Liver")
# Data available for compounds in Rat in-vitro liver Repeat data
tggates_compounds(species = "Rat", data_type = "in_vitro", tissue = "Liver", dose_type = "Repeat")
# Data available for compounds in human in-vitro liver Single data
tggates_compounds(species = "Human", data_type = "in_vitro", tissue = "Liver")
# Data available for compounds in Rat in-vivo kidney single data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "kidney", dose_type = "single")
# Data available for compounds in Rat in-vivo kidney Repeat data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "kidney", dose_type = "Repeat")
}
\seealso{
Other Open TG-GATEs Data Download Helpers: 
\code{\link{tggates_link}()}
}
\concept{Open TG-GATEs Data Download Helpers}
