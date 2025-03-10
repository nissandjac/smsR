
<!-- README.md is generated from README.Rmd. Please edit that file -->

# smsR

<!-- badges: start -->
<!-- badges: end -->

smsR is a seasonal stock assessment model package using TMB to estimate
parameters. The model is the current assessment model for sandeel (4
stocks) and sprat in the North Sea.

## Installation

You can install the development version of smsR as

``` r
# Anonymized for review 
```

## Using the model

The model can be run with the provided data set `sandeel_1r` and the
Sandeel 1r standard configuration as

``` r
library(smsR)


# Set up sandeel for area 1r

df.tmb <- get_TMB_parameters(
  mtrx = sandeel_1r$lhs, # List that contains M, mat, west, weca
  Surveyobs = sandeel_1r$survey, # Survey observations
  Catchobs = sandeel_1r$Catch, # Catch observations
  years = 1983:2021, # Years to run
  nseason = 2, # Number of seasons
  useEffort = TRUE, # Use effort to calculate F
  ages = 0:4, # Ages of the species
  recseason = 2, # Season where recruitment occurs
  CminageSeason = c(1, 1), # Minimum catch age per season
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = c(0, 1), # Minimum age in surveys
  Qmaxage = c(1, 3), # Max age in surveys
  Fbarage = c(1, 2), # Age use to calculate Fbar
  effort = sandeel_1r$effort, # Effort input
  blocks = c(1983, 1999), # Blocks with unique selectivity
  nocatch = sandeel_1r$nocatch, # Seasons where F is not calculated
  surveyStart = c(0.75, 0), #
  surveyEnd = c(1, 0), #
  surveySeason = c(2, 1), #
  surveySD = list(c(0, 1), c(1, 2)),
  catchSD = list(c(1, 3), c(1, 3)), # Catch CV groupings
  estSD = c(0, 2, 0), # Estimate CVs for 1) survey, 2) catch, 3) Stock recruitment relationship
  beta = 105809, # Hockey stick break point
  nllfactor = c(1, 1, 0.05) # Factor for relative strength of log-likelihood
)


parms <- getParms(df.tmb)

sas <- runAssessment(df.tmb, parms)
#> Congrats, your model converged
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

The resulting fit to the data can then be visualized in a summary output
showing SSB, Catch, recruitment and Fbar as

``` r
plot(sas, printFig = FALSE)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Or to evaluate the retrospective patterns as

``` r
mr <- mohns_rho(df.tmb, peels = 5, parms = parms, plotfigure = FALSE)

mr$p1()
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
