# Run a the 2022 stock assessment of sandeel in Area 1a
library(smsR)

years <- 1983:2021
nseason <- 2
ages <- 0:4

Catchobs <- df_to_matrix(sandeel_1r[['canum']], season =1:nseason)
Surveyobs <- survey_to_matrix(sandeel_1r[['survey']], year = years)
Feffort <- sandeel_1r[['effort']]
nocatch <- sandeel_1r[['nocatch']]

df.tmb <- get_TMB_parameters(
  mtrx = sandeel_1r[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  nseason = nseason, # Number of seasons
  useEffort = 1,
  ages = ages, # Ages of the species
  Fbarage = 1:2, # min and max age for Fbar calculations
  recseason = 2, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = c(0,1), # minimum age in surveys
  Qmaxage = c(1,3),
  isFseason = c(1,0), # Seasons to calculate fishing in
  effort = Feffort,
  powers =list(NA, NA),
  blocks = c(1983, 1989, 1999, 2005, 2010), # Blocks with unique selectivity
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = c(0, .75),
  surveyEnd =  c(0,1), # Does the survey last throughout the season it's conducted?
  surveySeason = c(2,1), # Which seasons do the surveys occur in
  surveyCV =  list(c(0,1), # Add ages for survey CV
                   c(1,2,3)),
  catchCV = list(c(0,1,3),
                 c(0,1,3)),
  recmodel = 2, # Chose recruitment model (2 = estimated)
  estCV = c(0,2,0), # Estimate
 # betaSR = 110000, # Hockey stick breakpoint (x-axis)
  nllfactor = c(1,1,1) # Factor for relative strength of log-likelihood

)



# Get initial parameter structure
parms <- getParms(df.tmb)


# Run the assessment

sas <- runAssessment(df.tmb, parms = parms)

# Plot standard graphs
p <- smsPlots(df.tmb, sas)

# See Mohns rho
m <- mohns_rho(df.tmb,parms, peels = 5, plotfigure = TRUE)
print(m$mohns)
# Show some other functions
require(ggplot2)

# Get SSB
SSB <- getSSB(df.tmb, sas)

# Recruitment
R <- getR(df.tmb, sas)

# N
N <- getN(df.tmb, sas)

ggplot(N, aes(x = years, y = N, color = factor(season), fill = factor(season)))+
  geom_line()+facet_wrap(~ages, scales = 'free_y')+
  geom_ribbon(aes(ymin = minSE, ymax = maxSE), alpha = .2, linetype = 0)+
  theme_classic()
#

# Catch
C <- getCatch(df.tmb, sas)

# Catch numbers
CN <- getCatchN(df.tmb, sas)

# Get all derived variables
sms.summary <- getSummary(df.tmb, sas)



# Get Fishing mortality
F0 <- getF(df.tmb, sas)
Fbar <- getFbar(df.tmb, sas)


# Get the CVs
surveyCV <- getSurveyCV(sas)

catchCV <- getCatchCV(sas)


# Residuals

residualCatch <- getResidCatch(df.tmb, sas)

ggplot(residualCatch, aes(x = years, y = ResidCatch, color = factor(ages)))+
  geom_col()+facet_grid(season~ages, scales = 'free_y')+
  theme_classic()+theme(legend.position = 'none')+
  geom_hline(aes(yintercept = 0), linetype = 2)

residualSurvey <- getResidSurvey(df.tmb, sas)

ggplot(residualSurvey, aes(x = years, y = ResidSurvey, color = factor(ages)))+
  geom_col()+facet_grid(survey~ages, scales = 'free_y')+
  theme_classic()+theme(legend.position = 'none')+
  geom_hline(aes(yintercept = 0), linetype = 2)


# Print the AIC of the model
AIC.sms(sas)

# Print the negative log likelihood
sas$opt[['objective']]

# Get estimated parameters
p.est <- getEstimatedParms(sas)
p.est

