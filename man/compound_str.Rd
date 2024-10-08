% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-structure&summary-matrix.R
\name{compound_str}
\alias{compound_str}
\title{Structure of levels in the data for a compound}
\usage{
compound_str(compound, metadata)
}
\arguments{
\item{compound}{Compound name (or abbreviation for \code{TG-GATEs} database)}

\item{metadata}{Metadata for given compounds, usually in the form of a data frame or
tibble object.}
}
\value{
A list of dose, time and replication
}
\description{
The \code{compound_str()} function is used to find the total number of dose levels,
the time points within each dose, and the number of replicates used for each
dose-time combination.

The \code{compound_str()} function calculates the total number of dose levels used for a
compound in the experiment, as well as the total number of time points for each dose
and the total number of replicates for each time point within each dose. For instance,
if a compound has 3 dose labels, there are 3 labels for the dose. If each dose label
is measured at 4 time points, then the number of labels for time levels would be
(4, 4, 4) corresponding to the 3 dose labels. Finally, if each dose-time combination
uses 3 replicates, then the number of labels for replicates would be 3 * 4 = 12
(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3). Therefore, the total number of labels at each
level is a product of the labels from earlier levels.
}
\examples{
sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
compound_str(compound = "Compound1",  metadata = sim_data$metadata)
}
