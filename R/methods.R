##Experiment
#' @export
'[[.Experiment' <- function(x, i, name = TRUE, drop = TRUE) x$get(i)
#' @export
detector_names.Experiment <- function(x) return(x$detector_names)
#' @export
plate_names.Experiment <- function(x) return(x$plate_names)
#' @export
sample_names.Experiment <- function(x) return(x$sample_names)
#' @export
n_samples.Experiment <- function(x) length(sample_names(x))
#' @export
n_plates.Experiment <- function(x) length(plate_names(x))
#' @export
n_detectors.Experiment <- function(x) length(detector_names(x))
#' @export
'[.Experiment' <- function(x, detectors, samples, plates, name = TRUE, drop = TRUE){
    # TODO: Warn when detectors/samples/plates not found.
    s <- Sys.time()
    plate_list <- x[['plates']]
    detector_idx = c(T)
    sample_idx = c(T)
    plate_idx = c(T)
    # Filter plates first...
    if(!missing(plates)){
        pn <- plate_names(x)
        stopifnot(all(plates %in% pn))
        plate_idx <- pn %in% plates
    }
    plate_list <- plate_list[plate_idx]
    new_plate_list <- list()
    for(plate_idx in seq_along(plate_list)){
        mod_plate <- plate_list[[plate_idx]]
        details <- list('BarCode' = plate_name(mod_plate))
        if(!missing(detectors)){
            detector_idx <- detector_names(mod_plate) %in% detectors
        }
        if(!missing(samples)){
            sample_idx <- sample_names(mod_plate) %in% samples
        }
        well_list <- mod_plate[['wells']][detector_idx & sample_idx]
        new_plate_list[[names(plate_list)[plate_idx]]] <- Plate$new(well_list, details=details)
    }
    newx = Experiment$new(plates=new_plate_list, input_directory=x[['input_directory']], history=x[['history']])
    updateHistory(newx, timepoint0=s, message=glue::glue("Subset to {n_plates(newx)} plates, {n_samples(newx)} samples, {n_detectors(newx)} detectors"))
    return(newx)
}

#' @export
print.Experiment <- function(x, ...){
    cat("EDS Experiment with: \n")
    cat("\t", glue::glue("{n_plates(x)} Plates"), "\n")
    cat("\t", glue::glue("{n_samples(x)} Samples"), "\n")
    cat("\t", glue::glue("{n_detectors(x)} Detectors"), "\n")
    cat('History: \n')
    print(x$get('history'))
    return(invisible(x))
}


## Plate
#' @export
print.Plate <- function(x, ...){
    cat("EDS Plate with: \n")
    cat("\t", glue::glue("{length(x[['wells']])} Wells"), "\n")
    cat("\t", glue::glue("{n_samples(x)} Samples"), "\n")
    cat("\t", glue::glue("{n_detectors(x)} Detectors"), "\n")
    cat("Metadata: \n")
    cat("Plate Name:", plate_name(x), "\n")
    cat("Input File:", plate_path(x), "\n")
}

#' @export
'[[.Plate' <- function(x, i, name = TRUE, drop = TRUE) x$get(i)
#' @export
sample_names.Plate <- function(x) return(x$sample_names)
#' @export
n_samples.Plate <- function(x) return(length(unique(sample_names(x))))
#' @export
plate_name.Plate <- function(x) return(x$plate_name)
#' @export
plate_path.Plate <- function(x) return(x$plate_path)
#' @export
detector_names.Plate <- function(x) return(x$detector_names)
#' @export
n_detectors.Plate <- function(x) return(length(unique(detector_names(x))))

## Well

#' @export
'[[.Well' <- function(x, i, name = TRUE, drop = TRUE) x$get(i)
#' @export
detector_name.Well <- function(x) return(x[['Detector']])
#' @export
sample_name.Well <- function(x) return(x[['Sample']])
#' @export
plate_name.Well <- function(x) return(x[['Plate']])



## Generics
#' @export
sample_names <- function(x) UseMethod("sample_names", x)
#' @export
detector_names <- function(x) UseMethod("detector_names", x)
#' @export
detector_name <- function(x) UseMethod("detector_name", x)
#' @export
sample_name <- function(x) UseMethod("sample_name", x)
#' @export
plate_name <- function(x) UseMethod("plate_name", x)
#' @export
plate_names <- function(x) UseMethod("plate_names", x)
#' @export
n_samples <- function(x) UseMethod("n_samples", x)
#' @export
n_detectors <- function(x) UseMethod("n_detectors", x)
#' @export
n_plates <- function(x) UseMethod("n_plates", x)
#' @export 
plate_path <- function(x) UseMethod("plate_path", x)


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
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @export
updateHistory <- function(x, timepoint0, message) UseMethod("updateHistory", x)
#' @export
updateHistory.Experiment <- function(x, timepoint0, message){
    x$updateHistory(timepoint0, message)
    return(x)
}
