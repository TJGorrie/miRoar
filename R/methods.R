#' Update History of miRoar class object
#'
#' Update History of miRoar class object
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

updateHistory.miRoar <- function(x, timepoint0, message){
    x$history <- rbind(x$history, 
                        data.frame('Submitted' = timepoint0, 
                        'Finished' = Sys.time(), 
                        'Comment' = message))
    return(x)
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
#' @export
deltaCt.miRoar <- function(ct, 
                    method = c('global', 'endogenous','geNorm', 'NormFinder'), 
                    HKs = NULL, 
                    group = NULL,
                    ...
                    ){
    t0 <- Sys.time()
    # These should be done outside!
    ct2 <- ct[['CtAvg']]
    amp <- ct[['CrtAmp']]
    is.na(ct) <- ct2 <= 0 | ct2 >= 40
    is.na(ct) <- amp == -1
    ct3 <- ct2[!duplicated(rownames(ct2)),]
    dct <- deltaCt(ct3, method = method, HKs = HKs, group = group, ...)

    ct[['dCT']] <- dct
    ct <- updateHistory(ct, t0, sprintf('Calculate deltaCT using %s', method))
    return(ct)
}
