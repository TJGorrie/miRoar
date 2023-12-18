input_path <- "/Users/tyler/Downloads/Schizophrenia_miRNA_data/RE-ANALYSED"

#' @export
generateHistoryLine <- function(timepoint0, message){
    data.frame('Submitted' = timepoint0, 'Finished' = Sys.time(), 'Comment' = message)
}

#' @export
Experiment <- R6::R6Class("Experiment", 
    portable = TRUE,
    public = list(
        initialize = function(plates, input_directory, history=NULL){
            private[['plates']] <- plates 
            private[['input_directory']] <- input_directory
            if(!is.null(history)){
                private[['history']] <- history
            }
        },
        get = function(x, deep=FALSE) {
            if(deep) sapply(private[['plates']], function(plate) plate$get(x, deep=TRUE))
            else private[[x]]
        },
        updateHistory = function(timepoint0, message) {
            private$history <- as.data.frame(
                rbind(
                    private$history, 
                    data.frame(
                        'Submitted' = timepoint0, 
                        'Finished' = Sys.time(), 
                        'Comment' = message
                    )
                ),
                stringsAsFactors = FALSE
            )
        }
    ),
    private = list(
        "plates" = NA,
        "history" = NULL,
        "input_directory" = NA
    ),
    active = list(
        sample_names = function(value){
            if(missing(value)){
                unique(unlist(lapply(private$plates, sample_names)))
            } else {
                stop('Not Allowed!')
            }
        },
        detector_names = function(value){
            if(missing(value)){
                unique(unlist(lapply(private$plates, detector_names)))
            } else {
                stop('Not Allowed!')
            }
        },
        plate_names = function(value){
            if(missing(value)){
                unlist(lapply(private$plates, plate_name))
            } else {
                stop('Not Allowed!')
            }
        }
    )
)

#' @export
Plate <- R6::R6Class("Plate",
    portable = TRUE,
    public = list(
        initialize = function(x, details){
            private[['wells']] <- x
            private[['.plate_name']] <- details$BarCode
        },
        get = function(x, deep=FALSE) {
            if(deep){
                wd <- private[['wells']]
                sample_name <- sapply(wd, function(well) well$get('Sample'))
                detector_name <- sapply(wd, function(well) well$get('Detector'))
                well_values <- sapply(wd, function(well) well$get(x))
                if(length(well_values > length(sample_name))){
                    browser()
                    # Convert to array?
                } else {
                    # Must be an easier way, check SO.
                    df <- data.frame(x=well_values, 'sample_name' = sample_name, 'detector' = detector_name)
                    sdf <- split(df, df$sample_name)
                    data <- do.call('cbind', lapply(sdf, function(x) {y <- x$x; names(y) <- x$detector; return(y)}))
                }
                return(data)

            } else return(private[[x]])
        }
    ),
    private = list(
        "wells" = NA,
        ".plate_name" = NA
    ),
    active = list(
        sample_names = function(value){
            if(missing(value)){
                sapply(private$wells, sample_name)
            } else {
                stop("Stop")
            }
        }, 
        detector_names = function(value){
            if(missing(value)){
                sapply(private$wells, detector_name)
            } else {
                stop("Stop")
            }
        }, 
        plate_name = function(value){
            if(missing(value)){
                private$.plate_name
            } else {
                stop("Stop")
            }
        }
    )
)

#' @export
Well <- R6::R6Class("Well",
    # Philosophy. Wells should not be editted.
    portable=TRUE,
    public = list(
        initialize = function(x) {
            private[['Well Index']] <- x$Well
            private[['Detector']] <- x$Detector
            private[['Sample']] <- x$Sample
            private[['Raw Amp Score']] <- x$amp_score
            private[['OA Rox']] <- x$oa_rox
            private[['Task']] <- x$Task
            private[['Rn']] <- x$RnValues
            private[['deltaRn']] <- x$deltaRnValues
            # These need to be broken up...
            private[['cRT Matrix']] <- x$cRTValues
            # These should be explicitly named...
            for(i in names(x[['CTValues']])) private[[i]] <- x[['CTValues']][i]
        },
        get = function(x) return(private[[x]])
    ),
    private = list(
        "Well Index" = NA,
        "Detector" = NA,
        "Raw Amp Score" = NA,
        "Sample" = NA,
        "OA Rox" = NA,
        "Task" = NA,
        "Rn" = NA,
        "deltaRn" = NA,
        "cRT Matrix" = NA,
        "Ct" = NA,
        "Avg Ct" = NA,
        "Ct SD" = NA,
        "Delta Ct" = NA,
        "Qty" = NA,
        "Avg Qty" = NA,
        "Qty SD" = NA,
        "Amp Score" = NA,
        "Crt Amp" = NA,
        "Crt Raw" = NA,
        "Cq Conf" = NA,
        "Plate" = NA
    ),
    active = list(
        sample_name = function(value) if(missing(value)) private[['Sample']] else stop("Stop"),
        detector_name = function(value) if(missing(value)) private[['Detector']] else stop("Stop"),
        plate_name = function(value) if(missing(value)) private[['Plate']] else stop("Stop")
    )
)