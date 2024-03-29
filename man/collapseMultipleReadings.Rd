% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{collapseMultipleReadings}
\alias{collapseMultipleReadings}
\title{Collapse multiple of readings miRoar class object or Ct matrix}
\usage{
collapseMultipleReadings(x, method = c("mean", "median"), na.rm = T)
}
\arguments{
\item{x}{a Crt matrix or miRoar object}

\item{method}{Which method to compute average readings for}
}
\value{
x
}
\description{
Collapses multiple readings into a single reading, reported as either the mean or median value
}
\examples{
# readEDSfile()

}
\author{
Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
}
