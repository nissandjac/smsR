# Output utilities


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
getSSB <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  SSB <- data.frame(SSB = sdrep[rep.values == 'SSB',1])
  SSB$SE <- sdrep[rep.values == 'SSB',2]
  SSB$minSE <- SSB$SSB-2*SSB$SE
  SSB$maxSE <- SSB$SSB+2*SSB$SE
  SSB$years <- c(years,max(years)+1)

  return(SSB)
  }

#' Get a data frame of recruitment
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of recruitment. Last year is based on the Stock recruitment relationship. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
getR <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  R <- data.frame(R = sdrep[rep.values == 'Rsave',1])
  R$SE <- sdrep[rep.values == 'Rsave',2]
  R$minSE <- R$R-2*R$SE
  R$maxSE <- R$R+2*R$SE
  R$years <- c(years,max(years)+1)

  return(R)
}


#' Get numbers at age from fitted model
#'
#' @param df.tmb list of input parameters
#' @param sas fitted smsR model
#'
#' @return
#' data frame of Numbers at age. Last year is based on the expected survival from terminal year. minSE and maxSE is the 95% confidence intervals
#' @export
#'
#' @examples
getN <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  N <- data.frame(N = sdrep[rep.values == 'Nsave',1])
  N$SE <- sdrep[rep.values == 'Nsave',2]
  N$minSE <- N$N-2*N$SE
  N$maxSE <- N$N+2*N$SE
  N$ages <- df.tmb$age
  N$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears+1))
  N$years <- rep(c(years,max(years)+1), each = df.tmb$nage)

  return(N)
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

  tmp <- data.frame(F0 = sdrep[rep.values == 'F0',1])
  tmp$SE <- sdrep[rep.values == 'F0',2]
  tmp$minSE <- tmp$F0-2*tmp$SE
  tmp$maxSE <- tmp$F0+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nseason*df.tmb$nage)

  F0 <- tmp

  return(F0)
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
getResidCatch <- function(df.tmb, sas){

  reps <- sas$reps

  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  tmp <- data.frame(ResidCatch = sdrep[rep.values == 'resid_catch',1])
  tmp$ResidCatch[tmp$ResidCatch == -99] <- NA

  tmp$SE <- sdrep[rep.values == 'resid_catch',2]
  tmp$SE[is.na(tmp$ResidCatch)] <- NA
  tmp$minSE <- tmp$ResidCatch-2*tmp$SE
  tmp$maxSE <- tmp$ResidCatch+2*tmp$SE
  tmp$ages <- df.tmb$age
  tmp$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*(df.tmb$nyears))
  tmp$years <- rep(years, each = df.tmb$nseason*df.tmb$nage)

  Resid_Catch <- tmp

  return(Resid_Catch)
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
getEstimatedParms <- function(sas){

  parms <- data.frame(parameter = sas$opt$par,
                      se = sas$reps$sd)


  return(Resid_Catch)
}







