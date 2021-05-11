# miRoar

The miRoar R package provides a simple framework for analysing RT-qPCR data. Specifically designed for taq-man Open Arrays for 754 miRs the code should be flexible enough to handle any configuration of open array plates providing the data is in a .eds format. We plan to introduce support for text files in the future.

## Installation
The miRoar R package can be installed directly from github by using the devtools R package 
```
BiocManager::install('NormqPCR')
devtools::install_github('tjgorrie/miRoar')
```

## Features
* Read in EDS files
* Calculate delta CT values using multiple methods
* QC data
* View various plots

## Quick-guide
```
library(miRoar)

# Read in Data
raw <- batchReadEDS('pathtofolderwithEDSfiles', dir = TRUE)

# Perform light filtering according to manufacturers defaults
process <- setBadSignalsToNA(raw, maxCT = 40, minCT = 0, ampVal = 0, conf.val = .8) # set conf.val to 0 if you want to skip them...

# Remove deduplicated observations - this step should be replaced with a method to collapse the readings from multiple eds files 
dedup <- subset(process, which(!duplicated(getRownames(process))) )

# Calculate delta Ct values using global normalisation and the Crt values
normalised <- deltaCt(dedup, which = 'Crt', method='global')

# Optional Remove miR that do not appear in a % of samples (1 = 100% of samples)
noSig <- removeBadSignals(normalised, perc = 1) # remove miRs that do not express in any samples!

noSig

# Finally Extract dCt values for statistical analysis 
dCrt <- noSig$dCt
dim(dCrt) # miRs along rows, samples along columns


```

## To Do
* S4 methods
* ddCT analysis
* More QC
* Various plotting functions.
