# paperOptBalGPPS
Package and simulations related to "Optimally Balanced Gaussian Process Propensity Scores for Estimating Treatment Effects" by Brian Vegetabile, Dan L. Gillen, and Hal Stern

- This repository serves as a snapshot of the "gpbalancer" package as of 3/7/2019 to provide consistent and stable reproducibility of results within the paper.  

## Installation Instructions

1. Clone or download this repository
2. Using RStudio, open the file `paperOptBalGPPS.Rproj` file
3. Compile the package (on a Mac this can be done with "command+shift+b"
    - This may require the `devtools` package

## Replicating the Simulation Section

The simulations that were run for the paper are located within the ~Simulations/ folder.  Each simulation was run with `n_sim=1000`.

- `01-atesim.R` - Provides simulations for comparing methods for estimating the ATE
- `02-atesim-tablesforpaper.R` - Compiles the results from 01-atesim.R for building tables
- `03-attsim-revision.R` - Provides simulations for comparing methods for estimating the ATT
- `04-attsim-tablesforpaper.R` - Compiles the results from 03-attsim-mixtures.R for building tables
- `05-timingresults.R` - Compares timing results that were provided in the discussion section

### Results from simulations within paper

Within the simulation folder are two other folders which contain the simulation results that were used to construct the tables for the paper.

 - ~Simulations/ate-simresults/
     - `2019-02-08-nonparametric_odd-atesim-results.rds`
     - `2018-02-08-nonparametric_even-atesim-results.rds`
 - ~Simulations/att-simresults/
     - `2019-02-11-attsim-multivariate-results-revision.rds`

## Recreating the Application

To recreate the application section results the file `dw99_replication.R` is provided that utilizes out GP methodology.  A secondary analysis re-evaluates the DW99 data using the models defined in that paper, but utilizing IPW estimation.  These can be recreated using `reanalyzing_ds99.R`

- Data for the replication is within the data folder. 
- Note: Propensity score estimates that were used in this paper were too large to store on Github. 
