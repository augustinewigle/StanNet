% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertToBUGSnetRun.r
\name{convertToBUGSnetRun}
\alias{convertToBUGSnetRun}
\title{Convert a StanNetRun object to a BUGSnetRun object}
\usage{
convertToBUGSnetRun(StanNetobj)
}
\arguments{
\item{StanNetobj}{A StanNetRun object to be converted to a BUGSnetRun object}
}
\value{
A BUGSnetRun object
}
\description{
Takes a StanNetRun object and converts the samples object to something of the same form as the samples object from nma.run,
to be compatible with nma.fit and other plotting functions from BUGSnet
}
