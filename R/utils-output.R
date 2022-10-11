# Output utilities


#' Get a data frame of the spawning stock biomass
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of spawning biomass and 95th confidence intervals and standard error
#' @export
#'
#' @examples
#' SSB <- getSSB(df.tmb, sas)
#' print(SSB)
getSSB <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  SSB <- data.frame(SSB = sdrep[rep.values == 'logSSB',1])
  SSB$SE <- (sdrep[rep.values == 'logSSB',2])
  SSB$minSE <- SSB$SSB-2*SSB$SE
  SSB$maxSE <- SSB$SSB+2*SSB$SE
  SSB$years <- c(years,max(years)+1)

  SSB$SSB <- exp(SSB$SSB)
  SSB$minSE <- exp(SSB$minSE)
  SSB$maxSE <- exp(SSB$maxSE)


  return(SSB)
}

#' Get a data frame of the spawning stock biomass
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of spawning biomass. minSE and maxSE is the 95% confidence intervals, se is the standard error
#' @export
#'
#' @examples
getBiomass <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  N <- array(sdrep[rep.values == 'logBiomass',1], dim = c(df.tmb$nage, df.tmb$nyears, df.tmb$nseason),
             dimnames = list(df.tmb$age, df.tmb$years, 1:df.tmb$nseason))
  Nse <- array(sdrep[rep.values == 'logBiomass',2], dim = c(df.tmb$nage, df.tmb$nyears, df.tmb$nseason),
               dimnames = list(df.tmb$age, df.tmb$years, 1:df.tmb$nseason))

  Biomass <- N
  BiomassSE <- as.data.frame.table(Nse)

  Biomass.df <- as.data.frame.table(Biomass)
  names(Biomass.df) <- c('ages','years','season','Biomass')
  Biomass.df$SE <- BiomassSE$Freq
  Biomass.df$minSE <- Biomass.df$Biomass-2*Biomass.df$SE
  Biomass.df$maxSE <- Biomass.df$Biomass+2*Biomass.df$SE
  Biomass.df$ages <- as.numeric(as.character(Biomass.df$ages))
  Biomass.df$years <- as.numeric(as.character(Biomass.df$years))

  Biomass.df <- Biomass.df %>% dplyr::select(Biomass, SE, minSE, maxSE, ages, season,years)
  Biomass.df <- Biomass.df[-which(Biomass.df$Biomass == 0),]
  Biomass.df$Biomass <- exp(Biomass.df$Biomass)
  Biomass.df$minSE <- exp(Biomass.df$minSE)
  Biomass.df$maxSE <- exp(Biomass.df$maxSE)

  return(Biomass.df)
}
#' Retrieve the total catch from a fitted smsR object
#'
#' @param df.tmb input parameters
#' @param sas fitted smsR object
#'
#' @return returns a data frame
#' @export
#'
#' @examples
#' catch <- getCatch(df.tmb, sas)
#' print(catch)
getCatch <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(Catch = sdrep[rep.values == 'logCatchtot',1])
  tmp$SE <- sdrep[rep.values == 'logCatchtot',2]
  tmp$minSE <- tmp$Catch-2*tmp$SE
  tmp$maxSE <- tmp$Catch+2*tmp$SE
  tmp$years <- years

  tmp$Catch <- exp(tmp$Catch)
  tmp$minSE <- exp(tmp$minSE)
  tmp$maxSE <- exp(tmp$maxSE)

  Catch <- tmp


  return(Catch)
}



#' Get a data frame of recruitment from a fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of recruitment. Last year is based on the Stock recruitment relationship. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
#' R <- getR(df.tmb, sas)
#' print(R)
getR <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  R <- data.frame(R = sdrep[rep.values == 'logRec',1])
  R$SE <- sdrep[rep.values == 'logRec',2]
  R$minSE <- R$R-2*R$SE
  R$maxSE <- R$R+2*R$SE
  R$years <- c(years,max(years)+1)


  R$R <- exp(R$R)
  R$minSE <- exp(R$minSE)
  R$maxSE <- exp(R$maxSE)



  return(R)
}


#' Get numbers at age from a fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of Numbers at age. Last year is based on the expected survival from terminal year. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
#' N <- getN(df.tmb, sas)
getN <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  N <- data.frame(N = sdrep[rep.values == 'logN',1])
  N$SE <- sdrep[rep.values == 'logN',2]
  N$minSE <- N$N-2*N$SE
  N$maxSE <- N$N+2*N$SE
  N$ages <- df.tmb$age
  N$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  N$years <- rep(years, each = df.tmb$nage)
  N <- N[-which(N$N == 0),]

  N$N <- exp(N$N)
  N$minSE <- exp(N$minSE)
  N$maxSE <- exp(N$maxSE)


  return(N)
}



#' Get the number of catch individuals per age
#'
#' @param df.tmb input parameters
#' @param sas fitted smsR model
#'
#' @return
#' @export
#'
#' @examples
#' CatchN <- getCatchN(df.tmb, sas)
#' print(CatchN)
getCatchN <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(CatchN = sdrep[rep.values == 'logCatchN',1])
  tmp$SE <- sdrep[rep.values == 'logCatchN',2]
  tmp$minSE <- tmp$CatchN-2*tmp$SE
  tmp$maxSE <- tmp$CatchN+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp <- tmp[-which(tmp$CatchN == 0),]


  tmp$CatchN <- exp(tmp$CatchN)
  tmp$minSE <- exp(tmp$minSE)
  tmp$maxSE <- exp(tmp$maxSE)

  return(tmp)
}





#' Get fishing mortality at age from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
getF <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(F0 = sdrep[rep.values == 'logF0',1])
  tmp$SE <- sdrep[rep.values == 'logF0',2]
  tmp$minSE <- tmp$F0-2*tmp$SE
  tmp$maxSE <- tmp$F0+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp$F0[tmp$F0 == 0] <- -Inf
  tmp$minSE[tmp$F0 == -Inf] <- -Inf
  tmp$maxSE[tmp$F0 == -Inf] <- -Inf

  tmp$F0 <- exp(tmp$F0)
  tmp$minSE <- exp(tmp$minSE)
  tmp$maxSE <- exp(tmp$maxSE)

  F0 <- tmp

  return(F0)
}

#' Title
#'
#' @param df.tmb list of input parameters
#' @param sas smsR fitted model
#' @param Fbarage ages to calculate Fbar
#'
#' @return
#' @export
#'
#' @examples
getFbar <- function(df.tmb, sas, Fbarage = NULL){

  if(is.null(Fbarage)){
    warning('provide ages to calculate Fbar')
    Fbarage <- df.tmb$age[df.tmb$age > min(df.tmb$CminageSeason)]
  }

  F0 <- getF(df.tmb, sas)
  tmp <-  F0[F0$ages %in% Fbarage,] %>%
    dplyr::group_by(years, ages) %>%
    dplyr::summarise(Fbar0 = sum(F0),
                     minSE0 = sum(minSE),
                     maxSE0 = sum(maxSE)) %>%
    dplyr::group_by(years) %>%
    dplyr::summarise(Fbar = mean(Fbar0),
                     minSE = mean(minSE0),
                     maxSE = mean(maxSE0)
    )


  return(tmp)
}

#' Get fishin mortality at age from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
#'
#' getResidCatch(df.tmb, sas)
#' print(getResidCatch)
#'
getResidCatch <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(ResidCatch = sdrep[rep.values == 'resid_catch',1])
  tmp$ResidCatch[tmp$ResidCatch == -99] <- NA

  # tmp$SE <- sdrep[rep.values == 'resid_catch',2]
  # tmp$SE[is.na(tmp$ResidCatch)] <- NA
  # tmp$minSE <- tmp$ResidCatch-2*tmp$SE
  # tmp$maxSE <- tmp$ResidCatch+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nseason*df.tmb$nage)

  tmp <- tmp[-which(is.na(tmp$ResidCatch)),]

  Resid_Catch <- tmp

  return(Resid_Catch)
}


#' Get survey residuals from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
getResidSurvey <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(ResidSurvey = sdrep[rep.values == 'resid_survey',1])
  tmp$ResidSurvey[tmp$ResidSurvey == -99] <- NA

  tmp$SE <- sdrep[rep.values == 'resid_survey',2]
  # tmp$SE[is.na(tmp$ResidSurvey)] <- NA
  # tmp$minSE <- tmp$ResidSurvey-2*tmp$SE
  # tmp$maxSE <- tmp$ResidSurvey+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp$survey <- rep(1:df.tmb$nsurvey, each = df.tmb$nage*(df.tmb$nyears))
  tmp <- tmp[-which(is.na(tmp$ResidSurvey)),]

  Resid_Survey <- tmp

  return(Resid_Survey)
}


#' Retrieve estimated parameters
#'
#' @param sas stock assessment output from smsR
#'
#' @return
#' data frame of estimated parameters and their SE
#' @export
#'
#' @examples
#'
getEstimatedParms <- function(sas){

  parms.estimated <- data.frame(value = sas$opt$par,
                      se = sqrt(diag(sas$reps$cov.fixed)),
                      parameter = names(sas$opt$par))


  return(parms.estimated)
}



#' Retrieve estimated parameters
#'
#' @param sas stock assessment output from smsR
#'
#' @return
#' data frame of estimated parameters and their SE
#' @export
#'
#' @examples
#'
getCatchCV <- function(sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(catchCV = sqrt(sdrep[rep.values == 'SD_catch2',1]))

  tmp$SE <- sdrep[rep.values == 'SD_catch2',2]
  tmp$minSE <- tmp$catchCV-2*tmp$SE
  tmp$maxSE <- tmp$catchCV+2*tmp$SE

  tmp$ages <- df.tmb$age

  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage)

  tmp <- tmp[-which(tmp$catchCV == 0),] # Remove the ones that are not caught


  CatchCV  = tmp
  return(tmp)
}


#' Retrieve estimated parameters
#'
#' @param sas stock assessment output from smsR
#'
#' @return
#' data frame of estimated parameters and their SE
#' @export
#'
#' @examples
#'
getSurveyCV <- function(sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(surveyCV = sdrep[rep.values == 'SDS',1])

  tmp$SE <- sdrep[rep.values == 'SDS',2]
  tmp$minSE <- tmp$surveyCV-2*tmp$SE
  tmp$maxSE <- tmp$surveyCV+2*tmp$SE

  tmp$ages <- df.tmb$age

  tmp$season <- rep(1:df.tmb$nsurvey, each = df.tmb$nage)

  tmp <- tmp[-which(tmp$surveyCV == 0),] # Remove the ones that are not caught

  surveyCV <- tmp
  return(surveyCV)
}

#' Calculate AIC from smsR stock assessment object
#'
#' @param sas Stock assesment object from smsR
#' @param p penalty on number of parameters, default = 2
#' @param n AICc sample size, default = Inf
#'
#' @return
#' @export
#'
#' @examples
AIC.sms <- function(sas, p=2, n=Inf){

    opt <- sas$opt

    k = length(opt[["par"]])
    nll = opt[["objective"]]
    aic.sms = p*k + 2*nll + 2*k*(k+1)/(n-k-1)
    return( aic.sms )
}


#' Get summary of derived variables from fitted smsR object
#'
#' @param df.tmb
#' @param sas
#' @param Fbarage
#'
#' @return
#' @export
#'
#' @examples
getSummary <-function(df.tmb, sas, Fbarage = NULL){

  SSB <- getSSB(df.tmb,sas)
  R <- getR(df.tmb, sas)
  Catch <- getCatch(df.tmb, sas)
  Fbar <- getFbar(df.tmb,sas, Fbarage)


  df.out <- data.frame(
    years = SSB$years,
    SSB = SSB$SSB,
    R = R$R,
    Catch = c(Catch$Catch,NA),
    Fbar = c(Fbar$Fbar,NA)
  )


return(df.out)
}



