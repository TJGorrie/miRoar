.safe.as.numeric <- function(x){
    # Set "" or " " to zero for convenience.
    x[(x == " " | x == "")] <- 0
    ifelse(is.na(as.numeric(x)), x, as.numeric(x))
}

.processAnalysisTable <- function(data, colheads){
    pieces <- strsplit(data, '\t')
    CTVals <- pieces[[1]]
    names(CTVals) <- colheads
    CTNums <- .safe.as.numeric(CTVals[-c(1,2,3,4)])
    names(CTNums) <- colheads[-c(1,2,3,4)]
    RnVals <- matrix(nrow=40, .safe.as.numeric(pieces[[2]][-1]))
    dRnVals <- matrix(nrow=40, .safe.as.numeric(pieces[[3]][-1]))
    # Need to give these names!!!!!!!
    cRTVals <- matrix(ncol=5, .safe.as.numeric(pieces[[4]][-1]))
    return(list(
        "Well" = CTVals[1],
        "Sample" = CTVals[2],
        "Detector" = CTVals[3],
        "Task" = CTVals[4],
        "CTValues" = CTNums,
        "RnValues" = RnVals,
        "deltaRnValues" = dRnVals, 
        "cRTValues" = cRTVals
    ))
}

#' @importFrom tools file_path_sans_ext
#' @importFrom glue glue
#' @importFrom XML xmlToList
.extractEDSFile <- function(in_path, out_path, verbose=FALSE){
    bn <- tools::file_path_sans_ext(basename(in_path))
    if(verbose) message(glue::glue("Extracting data from {in_path}"))
    unzip(in_path, exdir=file.path(out_path, bn))
    data_directory <- file.path(out_path, bn, 'apldbio', 'sds')
    #experiment_xmllist <- XML::xmlToList(file.path(data_directory, 'experiment.xml'))
    # Need to handle:
    #[1] "Name"                           "RunState"                      
    #[3] "FinalAnalysisCompleted"         "CreatedTime"                   
    #[5] "ModifiedTime"                   "RunStartTime"                  
    #[7] "RunEndTime"                     "PreReadState"                  
    #[9] "FileName"                       "Label"                         
    #[11] "Type" 
    # "CocktailWastePercentage"        "ChemistryType" "UMMFactor"                      "CommonSampleConcentration"     "CommonSampleConcentrationUnits" "TCProtocolMode"                 "DNATemplateType"                "InstrumentTypeId"              "BlockTypeID"                    "PlateTypeID"                   
    #"ExperimentProperty"             "ExperimentProperty"            
    # "ExperimentProperty"             "ExperimentProperty" 
    # Behavious gnarly because of NULLs
    #sample_df <- data.frame(t(sapply(experiment_xmllist[names(experiment_xmllist) == 'Samples'], '[', TRUE)))
    #detector_df <- data.frame(t(sapply(experiment_xmllist[names(experiment_xmllist) == 'Detectors'], '[', TRUE)))
    #filterdata_xmllist <- XML::xmlToList(file.path(data_directory, 'filterdata.xml'))[-1]
    # Contains 40 PlatePointData objects
        # Each Contains: "Stage"     "Cycle"     "Step"      "Point"     "PlateData"
        # Stage Appears to be an Int...
        # Cycle Appears to be an Int corresponding to cycle number
        # Step Appears to be an Int...
        # Point Appears to be Int...
        # PlateData Contains:
            # "Rows"      "Cols"      "WellData"  "Attribute" "Attribute" "Attribute" "Attribute" "Attribute" "Attribute" "Attribute"
            # Row and Col are obvious, issue is row/col posix not recorded?
            # WellData is \t separated values e.g. strsplit(filterdata_xmllist[[1]][['PlateData']][['WellData']], '\t')[[1]]
            # One of the attributes contains info about temperature.
            # filterdata_xmllist[[1]][['PlateData']][names(filterdata_xmllist[[1]][['PlateData']]) == 'Attribute']
    #analysis_protocol_xmllist <- XML::xmlToList(file.path(data_directory, 'analysis_protocol.xml'))
    # Contains things about the machine prescribed analysis e.g.
    # AmpThreshold used etc unlist(analysis_protocol_xmllist[[7]])
    amp_scores <- read.table(file.path(data_directory, 'amp_score.txt'), skip=2)
    # Check if all detectors in amp scores, and vice verse for integrity
    # WellIndex, DetectorName, AmpScore(0-?)
    #stopifnot(all(unique(amp_scores[,2]) %in% unlist(detector_df$Name)))
    # These make up the start of each well?
    oa_rox <- read.table(file.path(data_directory, 'oa_rox.txt'), skip=2)
    #stopifnot(all(amp_scores[,1] %in% oa_rox[,1]))
    plate_setup_xmllist <- XML::xmlToList(file.path(data_directory, 'plate_setup.xml'))
    #"Name"                "BarCode"             "Description"        
    #[4] "Rows"                "Columns"             "PlateKind"          
    #[7] "FeatureMap"          "FeatureMap"          "FeatureMap"         
    #[10] "FeatureMap"          "FeatureMap"          "FeatureMap"         
    #[13] "Wells"               "MultiZoneEnabled"    "LogicalZone"        
    #[16] "PassiveReferenceDye"
    # One featuremap, contains sample info, per well... important!

    # multicomponentdata.xml
    # Last file to check...

    # analysis_result.txt
    analysis_results <- readLines(file.path(data_directory, 'analysis_result.txt'))[-1]
    colheads <- strsplit(analysis_results[1], '\t')[[1]]
    well_starts <- which(grepl('^[0-9]+\t', analysis_results))
    well_data <- split(analysis_results[-1], rep(well_starts, each=4))
    processed_well_data <- lapply(well_data, .processAnalysisTable, colheads=colheads)
    names(processed_well_data) <- sapply(processed_well_data, '[', 'Well')
    for(i in seq_along(oa_rox[,1])){
        well_name <- as.character(oa_rox[i,1])
        if(well_name %in% names(processed_well_data)) processed_well_data[[well_name]][['oa_rox']] <- oa_rox[i,2]
    }
    for(i in seq_along(amp_scores[,1])){
        well_name <- as.character(amp_scores[i,1])
        if(well_name %in% names(processed_well_data)) processed_well_data[[well_name]][['amp_score']] <- amp_scores[i,3]
    }
    WellData <- sapply(processed_well_data, function(x) Well$new(x))
    plate <- Plate$new(WellData, plate_setup_xmllist)
    return(plate)
}

#' @importFrom glue glue
processEDSFiles <- function(input_path, verbose=TRUE){
    edsTarget <- file.path(input_path, '*.eds')
    allEds <- Sys.glob(edsTarget)[1:2]
    n_eds <- length(allEds)
    if(verbose) message(glue::glue("{n_eds} files found!"))
    ext_folder <- file.path(input_path, 'ext')
    if(file.exists(ext_folder)) stop(glue::glue("{ext_folder} already exists, will not proceed"))
    dir.create(ext_folder)
    Plates <- sapply(allEds, .extractEDSFile, out_path=ext_folder, verbose=verbose)
    unlink(ext_folder, recursive=TRUE)
    return(Plates)
}

#' @export
readEDSExperiment <- function(input_directory){
    s <- Sys.time()
    plates <- processEDSFiles(input_directory)
    out <- Experiment$new(plates, input_path)
    updateHistory(out, s, glue::glue("Create EDSExperiment with {n_plates(out)} plates, {n_samples(out)} samples, {n_detectors(out)} detectors"))  
    return(out)             
}