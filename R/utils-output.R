# Output utilities


#' Get a data frame of the spawning stock biomass
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of spawning stock biomass and 95\% confidence intervals.
#' SE is standard error of log SSB
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

  SSB <- data.frame(SSB = sdrep[rep.values == 'logSSB',1])
  SSB$SE <- (sdrep[rep.values == 'logSSB',2])
  SSB$low <- SSB$SSB-2*SSB$SE
  SSB$high <- SSB$SSB+2*SSB$SE
  SSB$years <- c(years,max(years)+1)

  SSB$SSB <- exp(SSB$SSB)
  SSB$low <- exp(SSB$low)
  SSB$high <- exp(SSB$high)


  return(SSB)
}

#' Get a data frame of the biomass of each age group
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame containing the biomass of each age group.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log Biomass
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
  Biomass.df$low <- Biomass.df$Biomass-2*Biomass.df$SE
  Biomass.df$high <- Biomass.df$Biomass+2*Biomass.df$SE
  Biomass.df$ages <- as.numeric(as.character(Biomass.df$ages))
  Biomass.df$years <- as.numeric(as.character(Biomass.df$years))

  Biomass.df <- Biomass.df %>% dplyr::select(Biomass, SE, low, high, ages, season,years)
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
#' @return
#' data frame containing the biomass of catch each year.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log Catch
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

  tmp <- data.frame(Catch = sdrep[rep.values == 'logCatchtot',1])
  tmp$SE <- sdrep[rep.values == 'logCatchtot',2]
  tmp$low <- tmp$Catch-2*tmp$SE
  tmp$high <- tmp$Catch+2*tmp$SE
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

  tmp <- data.frame(Q = sdrep[rep.values == 'logQsurv',1])
  tmp$Q[tmp$Q == 0] <- NA
  tmp <- exp(tmp)
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
#' data frame of recruitment.
#' Last year is based on the Stock recruitment relationship.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log R

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

  R <- data.frame(R = sdrep[rep.values == 'logRec',1])
  R$SE <- sdrep[rep.values == 'logRec',2]
  R$low <- R$R-2*R$SE
  R$high <- R$R+2*R$SE
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
#' data frame of Numbers at age.
#' Last year is based on the expected survival from terminal year.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log N

#' @export
#'
#' @examples
#' N <- getN(df.tmb, sas)
getN <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  N <- data.frame(N = sdrep[rep.values == 'logN',1])
  N$SE <- sdrep[rep.values == 'logN',2]
  N$low <- N$N-2*N$SE
  N$high <- N$N+2*N$SE
  N$ages <- df.tmb$age
  N$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  N$years <- rep(years, each = df.tmb$nage)

  if(min(N$N) == 0){
    N <- N[-which(N$N == 0),]
  }


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
#' data frame containing the numbers of individuals by age in the catch each year.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log N
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

  tmp <- data.frame(CatchN = sdrep[rep.values == 'logCatchN',1])
  tmp$SE <- sdrep[rep.values == 'logCatchN',2]
  tmp$low <- tmp$CatchN-2*tmp$SE
  tmp$high <- tmp$CatchN+2*tmp$SE
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
#' data frame of fishing mortality at age.
#' Last year is based on the expected survival from terminal year.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log F0

#' @export
#'
#' @examples
getF <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  tmp <- data.frame(F0 = sdrep[rep.values == 'logF0',1])
  tmp$SE <- sdrep[rep.values == 'logF0',2]
  tmp$low <- tmp$F0-2*tmp$SE
  tmp$high <- tmp$F0+2*tmp$SE
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
#' data frame containing the average fishing mortality each year.
#' low and high are the 95\% confidence intervals.
#' SE is standard error of log Fbar
#' @export
#'
#' @examples
getFbar <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)

  tmp <- data.frame(Fbar = sdrep[rep.values == 'logFavg',1])
  tmp$SE <- sdrep[rep.values == 'logFavg',2]
  tmp$low <- tmp$Fbar-2*tmp$SE
  tmp$high <- tmp$Fbar+2*tmp$SE
  tmp$years <- df.tmb$years
  tmp$Fbar[tmp$Fbar == 0] <- -Inf
  tmp$low[tmp$Fbar == -Inf] <- -Inf
  tmp$high[tmp$Fbar == -Inf] <- -Inf

  tmp$Fbar <- exp(tmp$Fbar)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)

  return(tmp)
}

#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return

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
#' data frame of fishing mortality at age. Last year is based on the expected survival from terminal year. low and high are the 95\% confidence intervals
#' @export
#'
#' @examples
getResidSurvey <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  tmp <- data.frame(ResidSurvey = sdrep[rep.values == 'resid_survey',1])
  tmp$ResidSurvey[tmp$ResidSurvey == -99] <- NA

  tmp$SE <- sdrep[rep.values == 'resid_survey',2]
  # tmp$SE[is.na(tmp$ResidSurvey)] <- NA
  # tmp$low <- tmp$ResidSurvey-2*tmp$SE
  # tmp$high <- tmp$ResidSurvey+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp$survey <- rep(dimnames(df.tmb$Surveyobs)$survey, each = df.tmb$nage*(df.tmb$nyears))
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
  rep.values <- rownames(sdrep)

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

  tmp <- data.frame(surveyCV = sdrep[rep.values == 'SDS',1])

  tmp$SE <- sdrep[rep.values == 'SDS',2]
  tmp$low <- tmp$surveyCV-2*tmp$SE
  tmp$high <- tmp$surveyCV+2*tmp$SE

  tmp$ages <- df.tmb$age

  tmp$survey <- rep(dimnames(df.tmb$Surveyobs)$survey, each = df.tmb$nage)

  #tmp <- tmp[-which(tmp$surveyCV == 0),] # Remove the ones that are not caught

  surveyCV <- tmp
  return(surveyCV)
}

#' Calculate AIC from smsR stock assessment object
#'
#' @param object Stock assesment object from smsR
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

#' Make sure no boundaries are overrun in the model
#'
#' @param df.tmb
#' @param parms
#' @param mps
#'
#' @return
#' @export
#'
#' @examples
checkBoundaries <- function(df.tmb, parms, mps){







}



#' Get the stock recruitment relationship
#'
#' @param df.tmb data frame with input data
#' @param sas fitted smsR model
#'
#' @return
#' @export
#'
#' @examples
getSR <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)

  years <- df.tmb$years
  SSB <- getSSB(df.tmb, sas)

  alpha <-  sdrep[rep.values == 'logalpha',1]
  SE <- (sdrep[rep.values == 'logalpha',2])
  low <- alpha-2*SE
  high <- alpha+2*SE
  SSBrange <- seq(1, max(SSB$SSB), length.out = 100)

  SR <- exp(alpha)+log(SSBrange)
  SR[SSBrange > df.tmb$betaSR] <- exp(alpha)+log(df.tmb$betaSR)

  mins <- exp(low)+log(SSBrange)
  mins[SSBrange > df.tmb$betaSR] <- exp(low)+log(df.tmb$betaSR)
  maxs <- exp(high)+log(SSBrange)
  maxs[SSBrange > df.tmb$betaSR] <- exp(high)+log(df.tmb$betaSR)


  SR.out <- data.frame(SR = exp(SR), minSR = exp(mins), maxSR = exp(maxs) ,SSB = SSBrange )

  return(SR.out)
}


#' Get the fisheries selectivity
#'
#' @param df.tmb
#' @param sas
#'
#' @return
#' @export
#'
#' @examples
getSelex <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  tmp <- data.frame(Fsel = sdrep[rep.values == 'Fsel',1])
  tmp$SE <- sdrep[rep.values == 'Fsel',2]
  tmp$low <- tmp$Fsel-2*tmp$SE
  tmp$high <- tmp$Fsel+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)
  tmp$Fsel[tmp$Fsel == 0] <- -Inf
  tmp$low[tmp$Fsel == -Inf] <- -Inf
  tmp$high[tmp$Fsel == -Inf] <- -Inf

  tmp$Fsel <- exp(tmp$Fsel)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)
  tmp$blocks <- NA

  # Add blocks for plotting
  nblocks <- rep(NA, length(unique(df.tmb$bidx)))

  for(i in 1:length(nblocks)){
    if(i == 1){
      y1 <- df.tmb$years[1]
    }

    if(i == length(nblocks)){
      y2 <- df.tmb$years[df.tmb$nyears]
    }else{
      y2 <- df.tmb$years[which(df.tmb$bidx == i)[1]]
    }

    nblocks[i] <- paste(y1,y2, sep = '-')
    tmp$blocks[tmp$years%in% y1:y2] <- nblocks[i]

    y1 <- y2

  }

return(tmp)
}


#' Get the estimated survey
#'
#' @param df.tmb
#' @param sas
#'
#' @return
#' @export
#'
#' @examples
getSurvey <- function(df.tmb, sas){

  # Get the estimated survey
  reps <- sas$reps
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  tmp <- data.frame(surveyest = sdrep[rep.values == 'Surveyout',1])
  tmp$SE <- sdrep[rep.values == 'Surveyout',2]
  tmp$low <- tmp$surveyest-2*tmp$SE
  tmp$high <- tmp$surveyest+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$survey <- rep(dimnames(df.tmb$Surveyobs)$survey, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nage)

  tmp <- tmp %>% dplyr::filter(surveyest > 0)

  tmp$surveyest <- exp(tmp$surveyest)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)


  return(tmp)
}


getSelectivity <- function(df.tmb, sas){

  # Get the estimated survey
  reps <- sas$reps
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years

  tmp <- data.frame(selec = sdrep[rep.values == 'log_exp_pattern',1])
  tmp$SE <- sdrep[rep.values == 'log_exp_pattern',2]
  tmp$low <- tmp$selec-2*tmp$SE
  tmp$high <- tmp$selec+2*tmp$SE
  tmp$ages <- df.tmb$age

  tmp$selec[tmp$selec < -100] <- -Inf

  tmp$selec <- exp(tmp$selec)
  tmp$low <- exp(tmp$low)
  tmp$high <- exp(tmp$high)


  return(tmp)
}








