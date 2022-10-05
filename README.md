# smsR
## The Stochastic MultiSpecies stock assessment model (SMS)

This is the R-package for the stochastic multispecies model. Extended description incoming.

To run smsR you need Rtools with a version matching your current R version. The library depends primarily on the 'TMB' package to estimate parameters, but has other dependencies to calculate internal functions.

To install smsR run

```
require(devtools)
# Might take a while
install_github('https://github.com/nissandjac/smsR', dependencies = TRUE)

```

See (examples) folder for an example using sandeel in Area 1r in the North Sea
