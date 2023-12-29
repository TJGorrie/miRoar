# qPCRExperiment

The qPCRExperiment R package provides a simple framework for analysing RT-qPCR data. Originally this package was designed to analyse Taq-Man Open Arrays however the code should be flexible enough to handle any configuration of files that is in .eds format. 


## Installation
The qPCRExperiment R package can be installed directly from github by using the devtools R package 
```
BiocManager::install('NormqPCR')
devtools::install_github('tjgorrie/miRoar')
```

## Current Features
* Read in EDS files
* Calculate delta CT values using multiple methods

# [Quick Guide](vignettes/qPCRExperiment.html)