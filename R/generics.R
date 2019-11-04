#' Update History of miRoar class object
#'
#' Update History generic
#'
#' @param x character, The file name to be read in
#' @param timepoint0 Logical, whether or not some message are output
#' @param message Logical, controlled by 
#' 
#' @return x
#'
#' @examples
#' # readEDS()
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @export
updateHistory <- function(x, timepoint0, message){
    UseMethod("updateHistory", x)
}