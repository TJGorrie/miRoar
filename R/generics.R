#' Collapse multiple of readings miRoar class object or Ct matrix
#'
#' Collapses multiple readings into a single reading, reported as either the mean or median value
#'
#' @param x a Crt matrix or miRoar object
#' @param method Which method to compute average readings for
#' @return x
#'
#' @examples
#' # readEDSfile()
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @export
collapseMultipleReadings <- function(x, method = c("mean", "median"), na.rm=T){
    UseMethod('collapseMultipleReadings', x)
}

#' Set Bad Signals to NA
#'
#' Set Bad Signals to NA
#'
#' @param x a miRoar object
#' @param maxCT a
#' @param minCT b
#' @param ampVal c
#' @param conf.val d
#' @return x
#'
#' @examples
#' # readEDSfile()
#'
#' @author Tyler Gorrie-Stone \email{tyler.gorrie-stone@diamond.ac.uk}
#' @export
setBadSignalsToNA <- function(x, maxCT = 40, minCT = 0, ampVal = 0, conf.val = .8){
    UseMethod('setBadSignalsToNA', x)
}
