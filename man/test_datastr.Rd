% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error-check.R
\name{test_datastr}
\alias{test_datastr}
\title{Validate Data Structure for ToxAssay}
\usage{
test_datastr(expr_str, error_call = rlang::caller_env())
}
\arguments{
\item{expr_str}{An object to be validated. It should be either an object of class \code{ToxAssay} or a list.}
}
\value{
If the input is valid, the function returns \code{TRUE}. If the input is invalid, an error message is triggered.
}
\description{
This function checks whether the provided data structure is a valid \code{ToxAssay} object or a list with the correct format.
It ensures that the input is either a \code{ToxAssay} object or a list of length 5, containing the elements \code{group}, \code{compound}, \code{dose}, \code{time}, and \code{replication}.
}
\examples{
# Example of a valid ToxAssay object or list
valid_list <- list(group = "A", compound = "Methimazole", dose = "High", time = "24h", replication = 3)
test_datastr(valid_list)

# Example of an invalid list (less than 5 elements)
invalid_list <- list(group = "A", compound = "Methimazole")
# This will trigger an error
\dontrun{
test_datastr(invalid_list)
}
}
