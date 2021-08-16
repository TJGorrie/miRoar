.readEDSFile <- function(file){
    bn <- gsub('.eds', '', file)
    message(sprintf('Reading %s.', bn))
    unzip(file, exdir=bn)
    newpath <- sprintf('%s/apldbio/sds/', bn)
    read_lines <- readLines(sprintf('%s/analysis_result.txt', newpath))
    tabbed_lines <- lapply(lapply(read_lines, strsplit, split='\t'), '[[', 1)
    headers <- tabbed_lines[c(1,2)][[2]]
    headers <- gsub(' ', '_', headers)
    headers <- paste0('.', headers)
    tabbed_lines <- tabbed_lines[-c(1,2)]
    #iterations <- seq(1, length(tabbed_lines), 3) # This is faster but not guaranteed to work...
    iterations <- which(sapply(tabbed_lines, length) == 15)
    pb <- progress::progress_bar$new(total = length(iterations))
    # amp_score.txt processing
    amp <- read.table(sprintf('%samp_score.txt', newpath), skip = 1,
                      stringsAsFactors = FALSE, sep = '\t', header = TRUE)
    rownames(amp) <- as.character(amp[,1])
    # oa_rox.txt
    rox <- read.table(sprintf('%soa_rox.txt', newpath), skip = 1,
                    stringsAsFactors = FALSE, sep='\t', header = TRUE)
    rownames(rox) <- as.character(rox[,1])

    # puredyematrix # Assume that well index is coordinate with an R matrix, i might be wrong but oh well...
    puredye <- readLines(sprintf('%spuredyematrix.txt', newpath))
    puredye <- do.call('rbind',lapply(puredye[-c(1:13, 46)], function(x) strsplit(x, split=',')[[1]])) 


    wells <- lapply(iterations, function(i){
        ctvals <- tabbed_lines[[i]]
        names(ctvals) <- headers
        Rnvals <- list(.Rn_values = as.numeric(tabbed_lines[[i+1]][-1]))
        dRnvals <- list(.Delta_Rn_values = as.numeric(tabbed_lines[[i+2]][-1]))
        cmv <- list(.Crt_Matrix_values = as.numeric(tabbed_lines[[i+3]][-1]))
        others <- list(
            .Amp.Score = as.numeric(amp[ctvals[1], 3]), 
            .ROX.Reading=as.numeric(rox[ctvals[1], 2]), 
            .Pure.Dye = as.numeric(puredye[as.numeric(ctvals[1])]), # I am assuming this is the correct index here... TS may have done things backwards!
            .filename = bn
        )
        values <- c(lapply(ctvals, asnumericorcharacter), Rnvals, dRnvals, cmv, others)
        pb$tick()
        Well$new(values)
    })

    # cleanup
    unlink(bn, recursive = TRUE, force = TRUE)

    return(wells)
}

readEDSFiles <- function(files){
    s <- Sys.time()
    data <- unlist(lapply(files, .readEDSFile))
    mirror <- miRoar$new(data, time.submitted=s)
    return(mirror)
}

asnumericorcharacter <- function(x){
    suppressWarnings(ifelse(any(is.na(as.numeric(x))), x, as.numeric(x)))
}