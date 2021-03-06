% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_LP.r
\name{compare_LP}
\alias{compare_LP}
\title{Compare the log-posterior of a model in R vs. Stan}
\usage{
compare_LP(stanObj, testPars, detail = F)
}
\arguments{
\item{stanObj}{A StanNetRun object created by running \code{runStan} with \code{run = TRUE}}

\item{testPars}{A list of named lists of parameters to use for testing}

\item{detail}{Logical, TRUE if a numerical value indicating result of the test should be returned (default), or FALSE if the values of the log-posteriors and their differences should be returned}
}
\value{
3 vectors containing value of the log posterior in R and Stan and their differences, OR value which indicates the outcome of the comparison:
0 if the log posteriors were equal up to the same constant (may include instances where both R and Stan evaluated the log posterior to be infinite)
-1 if the log posteriors were not equal up to the same constant
-2 if the parameter settings resulted in too many infinite values and the differences could not be compared to check for equality
}
\description{
Takes in a StanNetRun object and list of parameters and tests if the log-posteriors are equal (up to a constant) between the R and the Stan code
}
