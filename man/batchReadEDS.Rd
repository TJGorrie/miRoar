% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inout.R
\name{batchReadEDS}
\alias{batchReadEDS}
\title{Read single EDS file}
\usage{
batchReadEDS(EDSfiles, directory = FALSE, ...)
}
\arguments{
\item{EDSfiles}{character, The file name to be read in}

\item{directory}{Logical, if TRUE treat EDSfiles as directory and recursively looked for eds files within it}

\item{...}{Other arguements passed to dir}
}
\value{
The contents of many EDS files in a list containing:
}
\description{
Extracts and compiles a single EDS file
}
\examples{
# batchReadEDS()

}
\seealso{
\code{\link{readEDS}}
}
\author{
Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
}
\keyword{eds}
