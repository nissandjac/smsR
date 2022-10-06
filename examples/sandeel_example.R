# Run a the 2022 stock assessment of sandeel in Area 1a
library(smsR)

years <- 1983:2021
nseason <- 2
ages <- 0:4

names(sandeel_1r)

Catchobs <- df_to_matrix(sandeel_1r[['canum']], season =1:nseason)
Surveyobs <- survey_to_matrix(sandeel_1r[['survey']],year = years, season = 1:nseason)
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
  recseason = 2, # Season where recruitment occurs
  CminageSeason = c(0,0),
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = c(0,1), # minimum age in surveys
  Qmaxage = c(1,3),
  isFseason = c(1,0), # Seasons to calculate fishing in
  effort = Feffort,
  powers =list(NA),
  blocks = c(1983,1989, 1999,2005 ,2010), # Blocks with unique selectivity
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
  beta = 110000, # Hockey stick plateau
  nllfactor = c(1,1,0.05) # Factor for relative strength of log-likelihood

)


# Get initial parameter structure
parms <- getParms(df.tmb)


# Run the assessment

sas <- runAssessment(df.tmb, parms = parms)

p <- smsPlots(df.tmb, sas$reps, Fbarage = 1:2)


# Show some other functions
# Get SSB
SSB <- getSSB(df.tmb, sas)

# Recruitment
R <- getR(df.tmb, sas)

# N
N <- getN(df.tmb, sas)

ggplot2::ggplot(N[N$season == 1,], aes(x = years, y = N, color = factor(ages)))+
  geom_line()+facet_wrap(~ages, scales = 'free_y')+theme_classic()+theme(legend.position = 'none')

