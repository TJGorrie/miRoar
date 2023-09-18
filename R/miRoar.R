miRoar <- R6Class(
    'miRoar',
    private = list(
        .miRs = NA,
        .miR_well_indicies = NA,
        .samples = NA,
        .sample_well_indicies = NA,
        .history = NA
    ),
    public = list(
        Well_data = NA,
        initialize = function(x, time.submitted = Sys.time()){
            self$Well_data <- x

            detectors_vec <- sapply(x, '[[', 'Detector')
            unique_detectors <- unique(detectors_vec)
            names(unique_detectors) <- unique_detectors
            samples_vec <- sapply(x, '[[', 'Sample_Name')
            samples <- unique(samples_vec)
            names(samples) <- samples
            private$.miRs <-  unique_detectors
            private$.miR_well_indicies <- lapply(unique_detectors, function(mir) which(mir == detectors_vec))
            private$.samples <- samples
            private$.sample_well_indicies <- lapply(samples, function(sample) which(sample == samples_vec))

            private$.history <- data.frame('Submitted' = time.submitted, 'Finished' = Sys.time(), 'Comment' = 'Created miRoar object')
        },
        print = function(...){
            cat(sprintf("miRoar object consisting of %s miRs and %s Samples. \n", self$nrow, self$ncol))
            invisible(self)
        },
        nrow = function(...){
            length(private$.miRs)
        },
        ncol = function(...){
            length(private$.samples)
        },
        dim = function(...){
            c(self$nrow, self$ncol)
        },
        # Generate Matrix?
        get_data_matrix = function(
            type=c("Crt_Raw", "Ct", "Avg_Ct",
                "Pure.Dye", "ROX.Reading", "Amp.Score",
                "Cq_Conf",  "Crt_Amp", "Amp_Score", "Qty_SD", "Avg_Qty", 
                "Qty", "Delta_Ct", "Ct_SD","Task", "Well"
            ), 
            collapse_multiple_readings = c('median', 'mean', 'min', 'max', 'sd')
            ){
            type <- match.arg(type)
            collapse_multiple_readings <- match.arg(collapse_multiple_readings)
            collapse_multiple_readings <- switch(collapse_multiple_readings,
                'median' = median,
                'mean' = mean,
                'min' = min,
                'max' = max,
                'sd' = sd
            )
            # Get Tall Skinny Table
            data <- data.frame(
                Detector = sapply(self$Well_data, '[[', 'Detector'), 
                Sample_Name = sapply(self$Well_data, '[[', 'Sample_Name'), 
                Value = as.numeric(sapply(self$Well_data, '[[', type)),
                hide = sapply(self$Well_data, '[[', 'show')
            )
            data <- data[data$hide,]
            tibble <- tidyr::pivot_wider(data, names_from = Sample_Name, values_from = Value, values_fn = collapse_multiple_readings, values_fill=NA)
            rn <- tibble$Detector
            mat <- as.matrix(tibble[,-1])   
            rownames(mat) <- rn
            return(mat) 
        },
        get_data_array = function(...){
            # "Crt_Matrix_values", # "Delta_Rn_values", "Rn_values", 
            return(0)
        }
    ),
    active = list(
        miRs = function(value) if(missing(value)) private$.miRs else stop('`$miRs` is read only and cannot be changed!', call. = FALSE),
        samples = function(value) if(missing(value)) private$.samples else stop('`$samples` is read only and cannot be changed!', call. = FALSE),
        showing = function(value) if(missing(value)) sapply(self$Well_data, '[[', 'show'),
        history = function(value){
            if (missing(value)) private$.history 
            else {
                stopifnot(length(value) == 3)
                private$.history <- rbind(private$.history, value)    
            }
        },
        hidden_wells = function(value){
            vec <- sapply(self$Well_data, '[[', 'show')
            if (missing(value)) vec
            else {
                s = Sys.time()
                stopifnot(length(value) == length(vec))
                stopifnot(is.logical(value))
                for(i in seq_along(value)){
                    self$Well_data[[i]]$show <- value[i]
                }
            }
        }
    )
)

generateHistoryLine <- function(timepoint0, message){
    data.frame(Submitted = timepoint0, Finished = Sys.time(), Comment = message)
}
