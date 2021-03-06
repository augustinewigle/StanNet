% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rprepoislogpost.R
\name{rprepoislogpost}
\alias{rprepoislogpost}
\title{Calculate log-posterior in R}
\usage{
rprepoislogpost(testPars, stanObj)
}
\arguments{
\item{testPars}{A named list of parameters to use}

\item{stanObj}{A StanNetRun object to supply the data}
}
\value{
The value of the log-posterior distribution at the test parameters and data
}
\description{
Computes the log-posterior for a non-centred random effect model with poisson likelihood using the log link and default priors
}
