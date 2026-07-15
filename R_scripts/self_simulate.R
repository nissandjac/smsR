# Self simulate model #
library(smsR)
library(tidyverse)

# Get OM parameters #

pars <- sim_OM_parameters(nseason = 2,
                          Fpast = 0.2,
                          F0 = 0.2,
                          SDcatch = 0.2)

# Simulate OM

OM <- run.agebased.sms.op(pars)

# Fit the EM

df.tmb <- OM_to_smsR(OM,pars, recmodel = 2)

# Start with the true parameters
#df.tmb <- get_EM_parameters(pars, OM, df.tmb)
parms <- getParms(df.tmb)
# mps <- getMPS(tt, parms)
# Fit the estimation model
sas <- runAssessment(df.tmb, parms)

# Compare SSB
SSB <- getSSB(df.tmb, sas) %>% mutate(SSB_OM = c(OM$SSB, NA))

ggplot(SSB, aes(x = years, y = SSB))+geom_line()+
  geom_line(aes(y = SSB_OM), linetype = 2, col = 'red') +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray', alpha = .15)+
  theme_classic()


plotBubbles(sas)

plot(sas)


cobs <- getYield(df.tmb)
cest <- getCatch(df.tmb, sas)

plot(cobs)
lines(cest$Catch)
lines(OM$Catch, col = 'red')




