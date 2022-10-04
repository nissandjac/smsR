#' Title
#'
#' @param df.tmb list of tmb parameters 
#' @param reps estimated model. 
#' @param Fbarage Fbar age
#'
#' @return
#' @export
#'
#' @examples
#' smsplots(df.tmb, reps)
#' 
smsPlots <- function(df.tmb, reps, Fbarage = df.tmb$age[df.tmb$age>0]){
  
  
  # Plot SSB 
  
  # ages = 0:4
  # nseason = 2
  # 
  
  
  
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality 
  
  SSB <- data.frame(SSB = sdrep[rep.values == 'SSB',1])
  SSB$SE <- sdrep[rep.values == 'SSB',2]
  SSB$minSE <- SSB$SSB-2*SSB$SE
  SSB$maxSE <- SSB$SSB+2*SSB$SE
  SSB$model <- 'TMB'
  SSB$years <- c(years,max(years)+1)
  
  pssb <- ggplot(SSB, aes(x = years, y = SSB))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('SSB')+theme(legend.position = c(0.8,.8))
  
  # pssb
  # Recruitment 
  
  
  rec <- data.frame(rec = sdrep[rep.values == 'Rsave',1])
  rec$SE <- sdrep[rep.values == 'Rsave',2]
  rec$minSE <- rec$rec-2*rec$SE
  rec$maxSE <- rec$rec+2*rec$SE
  rec$model <- 'TMB'
  rec$years <- c(years,max(years)+1)
  
  # take care of crazy recruitment stuff 
  lims <- c(min(rec$rec)/2, max(rec$rec)*2)
  
  prec <- ggplot(rec, aes(x = years, y = rec))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('recruitment')+theme(legend.position = 'none')+coord_cartesian(ylim = lims)
  
  prec
  # Fishing mortality 
  
  quarters <- 1:nseason
  
  F0 <- data.frame(F0 = sdrep[rep.values == 'F0',1])
  F0$SE <- sdrep[rep.values == 'F0',2]
  F0$minSE <- F0$F0-2*F0$SE
  F0$maxSE <- F0$F0+2*F0$SE
  F0$age <- df.tmb$age
  F0$season <- rep(df.tmb$nseason, each = df.tmb$nage*df.tmb$nyears)
  F0$year <- rep(years, each = length(ages))

  Fbar <- F0[F0$age %in% Fbarage,] %>% group_by(year, age) %>% 
    summarise(Fbar0 = sum(F0),
              minSE0 = sum(minSE),
              maxSE0 = sum(maxSE)
    ) %>% group_by(year) %>% 
    summarise(Fmean = mean(Fbar0),
              minSE = mean(minSE0),
              maxSE = mean(maxSE0)
    ) %>% mutate(model = 'tmb')
  
  
  pF0 <- ggplot(Fbar, aes(x = year, y = Fmean))+geom_line(size = 1.3)+theme_classic()+
    geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('fishing mortality')+theme(legend.position = 'none')
  
  pF0
  
  
  # And catch 
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality 
  
  Catch <- data.frame(Catch = sdrep[rep.values == 'Catchtot',1])
  Catch$SE <- sdrep[rep.values == 'Catchtot',2]
  Catch$minSE <- Catch$Catch-2*Catch$SE
  Catch$maxSE <- Catch$Catch+2*Catch$SE
  Catch$model <- 'TMB'
  
  pCatch <- ggplot(Catch, aes(x = years, y = Catch))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('Catch')+theme(legend.position = c(0.8,.8))

  ls <- (pssb+pCatch)/(prec + pF0)
  print(ls)
  
  return(ls)
  }