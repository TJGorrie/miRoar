#' Well class
#' 
#' Describes a single Well and relevant properties
#' 
#' @import R6
#' @exportClass Well
#' @examples 
#' \dontrun{
#' example <- Well$new(name= '')
#' }
Well <- R6Class(
    'Well',
    private = list(
        .Well = NA,
        .Sample_Name = NA,
        .Detector = NA,
        .Task = NA,
        .Ct = NA,
        .Avg_Ct = NA,
        .Ct_SD = NA,
        .Delta_Ct = NA,
        .Qty = NA,
        .Avg_Qty = NA,
        .Qty_SD = NA,
        .Amp_Score = NA,
        .Crt_Amp = NA,
        .Crt_Raw = NA,
        .Cq_Conf = NA,
        .Rn_values = NA,
        .Delta_Rn_values = NA,
        .Crt_Matrix_values = NA,
        .Amp.Score = NA,
        .ROX.Reading = NA,
        .Pure.Dye = NA,
        .filename = NA,
        .show = TRUE,
        .norm_factor = 0
    ),
    public = list(
        initialize = function(values){
            sapply(names(values), function(x) private[[x]] <- values[[x]])
        },
        print = function(...){
            cat(sprintf("Sample: %s, miR: %s \n  CT: %s, ", self$Sample_Name, self$Detector, self$Avg_Ct))
            invisible(self)
        }
    ),
    active = list(
        Well = function(value) if(missing(value)) private$.Well else stop('`$Well` is read only and cannot be changed!', call. = FALSE),
        Sample_Name = function(value) if(missing(value)) private$.Sample_Name else stop('`$Sample_Name` is read only and cannot be changed!', call. = FALSE),
        Detector = function(value) if(missing(value)) private$.Detector else stop('`$Detector` is read only and cannot be changed!', call. = FALSE),
        Task = function(value) if(missing(value)) private$.Task else stop('`$Task` is read only and cannot be changed!', call. = FALSE),
        Ct = function(value) if(missing(value)) private$.Ct else stop('`$Ct` is read only and cannot be changed!', call. = FALSE),
        Avg_Ct = function(value) if(missing(value)) private$.Avg_Ct else stop('`$Avg_Ct` is read only and cannot be changed!', call. = FALSE),
        Ct_SD = function(value) if(missing(value)) private$.Ct_SD else stop('`$Ct_SD` is read only and cannot be changed!', call. = FALSE),
        Delta_Ct = function(value) if(missing(value)) private$.Delta_Ct else stop('`$Delta_Ct` is read only and cannot be changed!', call. = FALSE),
        Qty = function(value) if(missing(value)) private$.Qty else stop('`$Qty` is read only and cannot be changed!', call. = FALSE),
        Avg_Qty = function(value) if(missing(value)) private$.Avg_Qty else stop('`$Avg_Qty` is read only and cannot be changed!', call. = FALSE),
        Qty_SD = function(value) if(missing(value)) private$.Qty_SD else stop('`$Qty_SD` is read only and cannot be changed!', call. = FALSE),
        Amp_Score = function(value) if(missing(value)) private$.Amp.Score else stop('`$Amp_Score` is read only and cannot be changed!', call. = FALSE),
        Crt_Amp = function(value) if(missing(value)) private$.Crt_Amp else stop('`$Crt_Amp` is read only and cannot be changed!', call. = FALSE),
        Crt_Raw = function(value) if(missing(value)) private$.Crt_Raw else stop('`$Crt_Raw` is read only and cannot be changed!', call. = FALSE),
        Cq_Conf = function(value) if(missing(value)) private$.Cq_Conf else stop('`$Cq_Conf` is read only and cannot be changed!', call. = FALSE),
        Rn_values = function(value) if(missing(value)) private$.Rn_values else stop('`$Rn_values` is read only and cannot be changed!', call. = FALSE),
        Delta_Rn_values = function(value) if(missing(value)) private$.Delta_Rn_values else stop('`$Delta_Rn_values` is read only and cannot be changed!', call. = FALSE),
        Crt_Matrix_values = function(value) if(missing(value)) private$.Crt_Matrix_values else stop('`$Crt_Matrix_values` is read only and cannot be changed!', call. = FALSE),
        Amp.Score = function(value) if(missing(value)) private$.Amp.Score else stop('`$Amp.Score` is read only and cannot be changed!', call. = FALSE),
        ROX.Reading = function(value) if(missing(value)) private$.ROX.Reading else stop('`$ROX.Reading` is read only and cannot be changed!', call. = FALSE),
        Pure.Dye = function(value) if(missing(value)) private$.Pure.Dye else stop('`$Pure.Dye` is read only and cannot be changed!', call. = FALSE),
        show = function(value) if(missing(value)) private$.show else private$.show <- value,
        norm_factor = function(value) if(missing(value)) private$.norm_factor else private$.norm_factor <- value
    )
)



#mirror <- readEDSFiles(files)
#dat <- mirror$CT_Matrix(type='Ct', collapse_multiple_readings ='median')
#mirror$history
#mirror$history <- generateHistoryLine(Sys.time(), 'abc')


