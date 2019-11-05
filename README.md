# miRoar

The miRoar R package provides a simple framework for analysing RT-qPCR data. Specifically designed for taq-man Open Arrays for 754 miRs the code should be flexible enough to handle any configuration of open array plates providing the data is in a .eds format. We plan to introduce support for text files in the future.

## Installation
The miRoar R package can be installed directly from github by using the devtools R package 
```
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

raw <- readEDS('pathtofolderwithEDSfiles', dir = TRUE)



```

## To Do
* S4 methods
* ddCT analysis
* More QC