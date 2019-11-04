#' Read single EDS file
#'
#' Extracts and compiles a single EDS file
#'
#' @param file character, The file name to be read in
#' @param verbose Logical, whether or not some message are output
#' @param suppress Logical, controlled by 
#'
#' @return The contents of the EDS file in a list
#'
#' @examples
#' # readEDS()
#'

#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @seealso \code{\link{batchReadEDS}}
#' @keywords eds
#' @export
readEDS <- function(file, verbose = TRUE, suppress = FALSE){
    if(!suppress) message('Reading single EDS file, if you wish to use xy format please use batchReadEDS')
    
    # Test if file exists
    stopifnot(file.exists(file))
  
    # TODO: Test if file is an .eds file
    # Unzip File creating folder with name of file less .eds
    experiment <- gsub('.eds', '', file)
    if(verbose) message('Extracting: ', experiment)
    unzip(zipfile = file, exdir = experiment)
    # On Function exit - clean up files!
    # Structure is: apldbio/sds/
    newpath <- sprintf('%s/apldbio/sds/', experiment)

    # Processing Analysis.txt
    rl <- readLines(sprintf('%s/analysis_result.txt', newpath))
    # Detangle mixture of row lengths
    rlTab<- lapply(lapply(rl, strsplit, split='\t'), '[[', 1)
    rlTabSub12 <- rlTab[-c(1:2)]

    # Think of more elegant way to do this
    crtmat <- drn <- rn <- cts <- list()
    for(i in which(sapply(rlTabSub12, length) == 15)){
      cts[[length(cts)+1]] <- rlTabSub12[[i]]
      rn[[length(rn)+1]] <- rlTabSub12[[i+1]][-1]
      drn[[length(drn)+1]] <- rlTabSub12[[i+2]][-1]
      crtmat[[length(crtmat)+1]] <- rlTabSub12[[i+3]][-1]
    }
    cts <- as.data.frame(do.call('rbind', cts), stringsAsFactors = FALSE)
    colnames(cts) <- rlTab[[2]]
    rn <- do.call('rbind', rn)
    drn <- do.call('rbind', drn)
    crtmat <- do.call('rbind', crtmat)
    rownames(rn) <- rownames(drn) <- rownames(crtmat) <- cts$Well

    # Split data according to sample and wells?
    splitCts <- split(cts, cts[['Sample Name']])

    # amp_score.txt processing
    amp <- read.table(sprintf('%samp_score.txt', newpath), skip = 1,
                      stringsAsFactors = FALSE, sep = '\t', header = TRUE)
    rownames(amp) <- as.character(amp[,1])

    # oa_rox.txt
    rox <- read.table(sprintf('%soa_rox.txt', newpath), skip = 1,
                    stringsAsFactors = FALSE, sep='\t', header = TRUE)
    rownames(rox) <- as.character(rox[,1])

    # puredyematrix
    puredye <- readLines(sprintf('%spuredyematrix.txt', newpath))
    puredye <- do.call('rbind',lapply(puredye[-c(1:13, 46)], function(x) strsplit(x, split=',')[[1]]))


    splitCts <- lapply(splitCts, function(x){
        x$ampScore <- amp[x$Well,3]
        x$rox <- rox[x$Well,2]
        r <- rn[x$Well,]
        dr <- drn[x$Well,]
        crt <- crtmat[x$Well,]
        return(list(x, 'R' = r, 'dR' = dr, 'crtStats' = crt))
    })
    unlink(experiment, recursive = TRUE, force = TRUE)
    return(splitCts)
}

#' Read single EDS file
#'
#' Extracts and compiles a single EDS file
#'
#' @param EDSfiles character, The file name to be read in
#' @param directory Logical, if TRUE treat EDSfiles as directory and recursively looked for eds files within it
#' @param ... Other arguements passed to dir
#'
#' @return The contents of many EDS files in a list containing:
#'
#' @examples
#' # batchReadEDS()
#'
#' @author Tyler Gorrie-Stone \email{tgorri@essex.ac.uk}
#' @seealso \code{\link{readEDS}}
#' @keywords eds
#' @export
#' @exportClass miRoar
batchReadEDS <- function(EDSfiles, directory = FALSE, ...){
    s = Sys.time()
    # Read multiple EDS files
    # EDSfiles: a vector of eds files or a directory containing EDS files.
    if(directory) EDSfiles <- dir(EDSfiles, pattern = '.eds', ignore.case = TRUE, full.names = TRUE,...) # Allow for recursive through ...
    names(EDSfiles) <- gsub('.eds', '', basename(EDSfiles), ignore.case = TRUE)
    allData <- unlist(lapply(EDSfiles, readEDS, verbose = TRUE,
                           suppress = TRUE), recursive = FALSE)
    runname <- strsplit(names(allData), "[.]")
    diff.n <-  length(unique(sapply(lapply(allData, '[[', 1), nrow)))
    if(diff.n > 1){
        # perSample refactor
        allData <- lapply(allData, function(x){
            split <- strsplit(x[[1]][,'Detector'], '_')
            # Name handler...
            # 1 = 'BLANK' = remove!
            # 2 = 'manifest difference otherwise ######_mir - leave as is
            # 3 = 'Manifest difference: [species]-miR-###_######_mir, if 3 select last 2
            toKeep <- sapply(split, length) > 1
            newnames <- sapply(split, function(y){
                if(length(y) == 2) return(paste(y[1], y[2], sep='_'))
                else if(length(y) == 3) return(paste(y[2], y[3], sep='_'))
                return(y)
            })[toKeep]
            xclone <- lapply(x, function(z){
                return(z[toKeep,])
            })
            xclone[[1]][,'Detector'] <- newnames
            return(xclone)
        })
    } 

    get <- c('Ct', 'Avg Ct', 'Crt Amp', 'Crt Raw', 'rox', 'ampScore', 'Cq Conf')
    names(get) <- c('Ct', 'CtAvg', 'CrtAmp', 'Crt', 'rox', 'amp', 'Con') 

    biglist <- lapply(get, function(element, allData){
        sapply(allData, function(x){
            dat <- as.numeric(x[[1]][,element])
            names(dat) <- x[[1]][,'Detector']
            return(dat)
        })
    }, allData=allData)            

    # R and dR values in array Format
    biglist[['R']] <- array(as.numeric(unlist(lapply(allData, '[[', 'R'))),
                dim = c(nrow(allData[[1]][[1]]), ncol(allData[[1]][['R']]), length(allData)),
                dimnames = list(allData[[1]][[1]][,'Detector'], 1:ncol(allData[[1]][['R']]), names(allData))
                )
    biglist[['dR']] <- array(as.numeric(unlist(lapply(allData, '[[', 'dR'))),
                dim = c(nrow(allData[[1]][[1]]), ncol(allData[[1]][['dR']]), length(allData)),
                dimnames = list(allData[[1]][[1]][,'Detector'], 1:ncol(allData[[1]][['dR']]), names(allData))
                )
 
    #hkgmir <- names(which(table(rownames(biglist[[1]]))==16))
    biglist$FN <- sapply(runname, '[', 1)
    biglist$history <- data.frame('Submitted' = s, 'Finished' = Sys.time(), 'Comment' = 'Created object')
    class(biglist) <- 'miRoar'
    return(biglist)
}