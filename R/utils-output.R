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
  SSB$CV <- (sdrep[rep.values == 'logSSB',2])
  SSB$low <- SSB$SSB-2*SSB$CV
  SSB$high <- SSB$SSB+2*SSB$CV
  SSB$years <- c(years,max(years)+1)

  SSB$SSB <- exp(SSB$SSB)
  SSB$low <- exp(SSB$low)
  SSB$high <- exp(SSB$high)


  return(SSB)
}

#' Get a data frame of the spawning stock biomass
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of spawning biomass. low and high is the 95% confidence intervals, se is the standard error
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
  Biomass.df$CV <- BiomassSE$Freq
  Biomass.df$low <- Biomass.df$Biomass-2*Biomass.df$CV
  Biomass.df$high <- Biomass.df$Biomass+2*Biomass.df$CV
  Biomass.df$ages <- as.numeric(as.character(Biomass.df$ages))
  Biomass.df$years <- as.numeric(as.character(Biomass.df$years))

  Biomass.df <- Biomass.df %>% dplyr::select(Biomass, CV, low, high, ages, season,years)
  Biomass.df <- Biomass.df[-which(Biomass.df$Biomass == 0),]
  Biomass.df$Biomass <- exp(Biomass.df$Biomass)
  Biomass.df$low <- exp(Biomass.df$low)
  Biomass.df$high <- exp(Biomass.df$high)

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
  tmp$CV <- sdrep[rep.values == 'logCatchtot',2]
  tmp$low <- tmp$Catch-2*tmp$CV
  tmp$high <- tmp$Catch+2*tmp$CV
  tmp$years <- years

  tmp$Catch <- exp(tmp$Catch)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)

  Catch <- tmp


  return(Catch)
}


#' Get the survey catchability
#'
#' @param df.tmb sms input parameters
#' @param sas fitted model object
#'
#' @return
#' @export
#'
#' @examples
getCatchability <- function(df.tmb, sas){


  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(Q = sdrep[rep.values == 'Qsurv',1])
  tmp$age <- rep(df.tmb$age, df.tmb$nsurvey)
  tmp$survey <- rep(1:df.tmb$nsurvey, each = df.tmb$nage)



  return(tmp)
}

#' Retrieve the fishing selectivity
#'
#' @param df.tmb sms input parameters
#' @param sas fitted sms model
#'
#' @return
#' @export
#'
#' @examples
#' ## ' Not Run
#' sel <- getSel(df.tmb, sas)
#'
#' ## End
getSel<- function(df.tmb, sas){


  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  blocks <- unique(df.tmb$bidx)
  Fage <- df.tmb$Fminage:df.tmb$Fmaxage
  # Plot SSB, recruitment, catch and fishing mortality
  tmp <- data.frame(Fsel = 0, age = rep(df.tmb$age, max(df.tmb$bidx)+1),
                    block = rep(unique(df.tmb$bidx), each = df.tmb$nage))

  Fsel = array(as.numeric(exp(sdrep[rep.values == 'logFage',1])), dim = c(length(Fage), length(blocks)))

  for(i in 1:length(blocks)){
    tmp$Fsel[tmp$age %in% Fage & tmp$block == blocks[i]] <- Fsel[,i]/max(Fsel[,i])
    tmp$Fsel[tmp$age > max(Fage) & tmp$block == blocks[i]] <- Fsel[nrow(Fsel),i]/max(Fsel[,i])
  }


  return(tmp)
}


#' Get a data frame of recruitment from a fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of recruitment. Last year is based on the Stock recruitment relationship. low and high is the 95% confidence intervals
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
  R$CV <- sdrep[rep.values == 'logRec',2]
  R$low <- R$R-2*R$CV
  R$high <- R$R+2*R$CV
  R$years <- c(years,max(years)+1)


  R$R <- exp(R$R)
  R$low <- exp(R$low)
  R$high <- exp(R$high)



  return(R)
}


#' Get numbers at age from a fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of Numbers at age. Last year is based on the expected survival from terminal year. low and high is the 95% confidence intervals
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
  N$CV <- sdrep[rep.values == 'logN',2]
  N$low <- N$N-2*N$CV
  N$high <- N$N+2*N$CV
  N$ages <- df.tmb$age
  N$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  N$years <- rep(years, each = df.tmb$nage)
  N <- N[-which(N$N == 0),]

  N$N <- exp(N$N)
  N$low <- exp(N$low)
  N$high <- exp(N$high)


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
  tmp$CV <- sdrep[rep.values == 'logCatchN',2]
  tmp$low <- tmp$CatchN-2*tmp$CV
  tmp$high <- tmp$CatchN+2*tmp$CV
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp <- tmp[-which(tmp$CatchN == 0),]


  tmp$CatchN <- exp(tmp$CatchN)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)

  return(tmp)
}





#' Get fishing mortality at age from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. low and high is the 95% confidence intervals
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
  tmp$CV <- sdrep[rep.values == 'logF0',2]
  tmp$low <- tmp$F0-2*tmp$CV
  tmp$high <- tmp$F0+2*tmp$CV
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp$F0[tmp$F0 == 0] <- -Inf
  tmp$low[tmp$F0 == -Inf] <- -Inf
  tmp$high[tmp$F0 == -Inf] <- -Inf

  tmp$F0 <- exp(tmp$F0)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)

  F0 <- tmp

  return(F0)
}

#' Title
#'
#' @param df.tmb list of input parameters
#' @param sas smsR fitted model
#'
#' @return
#' @export
#'
#' @examples
getFbar <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(Fbar = sdrep[rep.values == 'logFavg',1])
  tmp$CV <- sdrep[rep.values == 'logFavg',2]
  tmp$low <- tmp$Fbar-2*tmp$CV
  tmp$high <- tmp$Fbar+2*tmp$CV
  tmp$years <- df.tmb$years
  tmp$Fbar[tmp$Fbar == 0] <- -Inf
  tmp$low[tmp$Fbar == -Inf] <- -Inf
  tmp$high[tmp$Fbar == -Inf] <- -Inf

  tmp$Fbar <- exp(tmp$Fbar)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)

  return(tmp)
}

#' Get fishin mortality at age from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. low and high is the 95% confidence intervals
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
  # tmp$low <- tmp$ResidCatch-2*tmp$SE
  # tmp$high <- tmp$ResidCatch+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)

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
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. low and high is the 95% confidence intervals
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
  # tmp$low <- tmp$ResidSurvey-2*tmp$SE
  # tmp$high <- tmp$ResidSurvey+2*tmp$SE
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



#' Retrieve estimated catch CV
#'
#' @param sas stock assessment output from smsR
#'
#' @return
#' data frame of estimated parameters and their SE
#' @export
#'
#' @examples
#'
getCatchCV <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(catchCV = sdrep[rep.values == 'SD_catch2',1])

  tmp$SE <- sdrep[rep.values == 'SD_catch2',2]
  tmp$low <- tmp$catchCV-2*tmp$SE
  tmp$high <- tmp$catchCV+2*tmp$SE

  tmp$ages <- df.tmb$age

  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage)

  #tmp <- tmp[-which(tmp$catchCV == 0),] # Remove the ones that are not caught


  CatchCV <- tmp
  return(CatchCV)
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
getSurveyCV <- function(df.tmb,sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(surveyCV = sdrep[rep.values == 'SDS',1])

  tmp$SE <- sdrep[rep.values == 'SDS',2]
  tmp$low <- tmp$surveyCV-2*tmp$SE
  tmp$high <- tmp$surveyCV+2*tmp$SE

  tmp$ages <- df.tmb$age

  tmp$survey <- rep(1:df.tmb$nsurvey, each = df.tmb$nage)

  #tmp <- tmp[-which(tmp$surveyCV == 0),] # Remove the ones that are not caught

  surveyCV <- tmp
  return(surveyCV)
}

#' Calculate AIC from smsR stock assessment object
#'
#' @param sas Stock assesment object from smsR
#' @param p penalty on number of parameters, default = 2
#' @param n AICc sample size, default = Inf
#' @aliases AIC AIC.sms
#'
#' @return
#' @export
#'
#' @examples
AIC.sms <- function(object, p=2, n=Inf, ...){

    opt <- object$opt

    k = length(opt[["par"]])
    nll = opt[["objective"]]
    aic.sms = p*k + 2*nll + 2*k*(k+1)/(n-k-1)
    return( aic.sms )
}

# Extract the log likelihood of a sms model
#
# @return object of class \code{logLik} with attributes
# \item{val}{log likelihood}
# \item{nobs,nall}{number of non NA observations initially supplied to TMB}
# \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.sms <- function(object, ...) {
  if(object$opt$convergence !=0){
    val <- NA
  }else val <- -object$opt$objective

  nobs <- nobs.sms(object)
  df <- length(object$opt[["par"]])
  structure(val, nobs = nobs, nall = nobs, df = df,
            class = "logLik")
}

#' @importFrom stats nobs
#' @export
nobs.sms <- function(object, ...) {
  sum(object$x$resid_survey !=-99) +sum(object$x$resid_catch !=-99)
}


#' Get summary of derived variables from fitted smsR object
#'
#' @param df.tmb
#' @param sas
#'
#' @return
#' @export
#'
#' @examples
getSummary <-function(df.tmb, sas){

  SSB <- getSSB(df.tmb,sas)
  R <- getR(df.tmb, sas)
  Catch <- getCatch(df.tmb, sas)
  Fbar <- getFbar(df.tmb, sas)


  df.out <- data.frame(
    years = SSB$years,
    SSB = SSB$SSB,
    minSSB = SSB$low,
    maxSSB = SSB$high,
    R = R$R,
    Catch = c(Catch$Catch,NA),
    Fbar = c(Fbar$Fbar,NA)
  )


return(df.out)
}




