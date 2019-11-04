#' Calculate delta Ct values from CT values using a variety of methods
#'
#' Calculate delta Ct values using global, endogenous, genorm or normfinder
#'
#' @param ct character, The file name to be read in
#' @param method Logical, whether or not some message are output
#' @param HKs Logical, controlled by 
#' @param group 
#' @param ...
#' 
#' @return deltaCT values
#'
#' @examples
#' # readEDS()
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @export
deltaCt <- function(ct, 
					method = c('global', 'endogenous','geNorm', 'NormFinder'), 
					HKs = NULL, 
					group = NULL,
					...
					){
   UseMethod("deltaCt", ct)
}

#' Calculate delta Ct values from CT values using a variety of methods
#'
#' Calculate delta Ct values using global, endogenous, genorm or normfinder
#'
#' @param ct character, The file name to be read in
#' @param method Logical, whether or not some message are output
#' @param HKs Logical, controlled by 
#' @param group 
#' @param ...
#' 
#' @return deltaCT values
#'
#' @examples
#' # readEDS()
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @seealso \code{\link{batchReadEDS}}
#' @keywords deltaCT CT genorm normfinder
#' @importFrom NormqPCR geomMean
#' @importFrom matrixStats colMedians
#' @export
deltaCt.matrix <- function(ct, 
					method = c('global', 'endogenous','geNorm', 'NormFinder'), 
					HKs = NULL, 
					group = NULL,
					...
					){
	method <- match.arg(method)

	sampleFactors <- switch(method, 
							'global' = colMedians(na.omit(ct)),
							'endogenous' = {
								stopifnot(is.character(HKs))
							    apply(as.matrix(ct[HKs, ]), 2, geomMean, na.rm = T)
							}, # handle for multiple HKs? geomMean?
							'geNorm' = {
								stopifnot(is.numeric(HKs))
								selectHKwrapper(ct, method = method, minNrHKs = HKs, group = group, ...)
								
							},
							'NormFinder' = {
								stopifnot(is.numeric(HKs))
								selectHKwrapper(ct, method = method, minNrHKs = HKs, group = group, ...)
								
							}
							)

	deltaCts <- sweep(ct, 2, sampleFactors)
	return(deltaCts)
}

#' A wrapper for the selectHKs function
#'
#' A wrapper to calculate the geometric mean from HK genes selected by genorm or normfinder algorithms
#'
#' @param x ct Values
#' @param method character, which method to use, either genorm or normfinder
#' @param minNrHKs Numeric, number of miRs
#' @param group factor, group variable for normFinder
#' @param trace IF messages are to be produced
#' @param log if values are in log space
#' @param na.rm if na values should be removed
#' 
#' @return the factor to be substract from CT values
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @seealso \code{\link{deltaCT}}
#' @keywords deltaCT CT genorm normfinder
#' @importFrom NormqPCR selectHKs geomMean
selectHKwrapper <- function(x, 
							method = c('geNorm', 'NormFinder'), 
							minNrHKs=2, group = NULL, trace=T, na.rm=T, log = T){
	genes <- selectHKs(t(x), group = group, method = method, minNrHKs = minNrHKs,
              	log = log, Symbols = rownames(x), trace = trace, na.rm = na.rm)
	res <- apply(x[genes[[1]][1:minNrHKs],], 2, geomMean, na.rm = TRUE)
	return(res)
}
