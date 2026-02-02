#
library(smsR)
year <- 1:30
nseason <- 2
# Required input
mtrx <- list(mat = sim_data$mat,
             M = sim_data$M,
             weca = sim_data$weca,
             west = sim_data$west
)
Surveyobs <- sim_data$Surveyobs
Catchobs <- sim_data$Catchobs
df.tmb <- get_TMB_parameters(mtrx,
                             Surveyobs, # From a simulated model
                             Catchobs, # From a simulated model
                             years = year,
                             nseason = nseason,
                             ages= 0:5,
                             recmodel = 3,
                             #randomR = 1
)
parms <- getParms(df.tmb)
mps <- getMPS(df.tmb, parms)
sas <- runAssessment(df.tmb, parms , mps = mps)
# Plot the fit
plot(sas)
