# paperOptBalGPPS
Package and simulations related to "Optimally Balanced Gaussian Process Propensity Scores for Estimating Treatment Effects" by Brian Vegetabile, Dan L. Gillen, and Hal Stern

- This repository serves as a snapshot of the "gpbalancer" package as of 4/12/2018 to provide consistent and stable reproducibility of results within the paper.  

## Installation Instructions

1. Clone or download this repository
2. Using RStudio, open the file paperOptBalGPPS.Rproj file
3. Compile the package (on a Mac this can be done with `command+shift+b'
    - This may require the devtools package

## Running Simulations

The simulations that were run for the paper are located within the ~Simulations/ folder

- 01-atesim.R - Provides simulations for comparing the ATE
- 02-atesim-tablesforpaper.R - Compiles the results from 01-atesim.R for building tables


