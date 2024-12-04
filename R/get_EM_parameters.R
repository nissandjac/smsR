#' Get parameters from an operating model
#'
#' @param df.OM operating model input
#' @param OM smsR fitted operating model
#' @param df.tmb data input for assessment model used to condition OM
#' @param obsSurvey should observed or modeled survey be used
#' @param obsCatch should observed or modeled catch be used
#'
#' @return returns a list with smsR data input
#' @export
#'
#' @examples get_EM_parameters(...)
get_EM_parameters <- function(df.OM,
                              OM,
                              df.tmb, # Fitted assessment model
                              obsSurvey = FALSE,
                              obsCatch = FALSE){



  nyear <- length(df.OM$years)

  if(obsSurvey == TRUE){
    # Surveyobs <- array(-1 , dim = c(df.assess$nage, nyear, df.assess$nsurvey))
    #
    # for(k in 1:df.assess$nsurvey){
    #   Surveyobs[,,k] <- df.assess$Surveyobs[,,k]
    # }
    Surveyobs <- df.OM$Surveyobs

  }else{
    Surveyobs <- OM$survey
  }

  Surveyobs[is.na(Surveyobs)] <- -1


  if(obsCatch == TRUE){
  Catchobs <- df.OM$Catchobs
  Catchobs[is.na(Catchobs)] <- -1
  }else{
  Catchobs <- apply(OM$CatchN.save.age, MARGIN = c(1,2,4), FUN = sum)
  }
  # Group index and effort #
  # Average over spatial domain

  weca <- apply(df.OM$weca, MARGIN = c(1,2,4), FUN = mean)
  west <- apply(df.OM$west, MARGIN = c(1,2,4), FUN = mean)
  mat <- apply(df.OM$mat, MARGIN = c(1,2,4), FUN = mean)
  M <- apply(df.OM$M, MARGIN = c(1,2,4), FUN = mean)

  mtrx <- list(weca = weca, west = west, mat = mat, M = M)


  df.sim <- get_TMB_parameters(
    mtrx = mtrx,
    Surveyobs = Surveyobs,
    Catchobs = Catchobs,
    propM = df.OM$propM,
    propF = df.OM$propF,
    years = df.OM$years,
    nseason = df.OM$nseason,
    nsurvey = nrow(Surveyobs),
    ages = df.OM$age,
    Fbarage = df.OM$Fbarage,
    recseason = df.OM$rseason,
    Fminage = df.tmb$Fminage,
    Fmaxage = df.tmb$Fmaxage,
    Qminage = df.tmb$Qminage,
    Qmaxage = df.tmb$Qmaxage,
    Qlastage = df.tmb$Qlastage,
    isFseason = df.tmb$isFseason,
    CminageSeason = df.tmb$CminageSeason,
    endFseason = df.tmb$endFseason,
    nocatch = df.tmb$nocatch,
    useEffort = df.tmb$useEffort,
    effort = df.tmb$effort,
    blocks = df.tmb$blocks,
    surveyStart = df.OM$surveyStart,
    surveyEnd = df.OM$surveyEnd,
    surveySeason =  df.OM$surveySeason,
    Mprior = df.tmb$Mprior,
    power = df.tmb$powers,
    Pred_in = df.tmb$Pred_in,
    surveyCV = df.tmb$surveyCV,
    catchCV = df.tmb$catchCVin,
    #MCV = df.tmb$Midx_CV Fix later
    estCV = df.tmb$estCV,
    betaSR = df.tmb$betaSR,
    nllfactor = df.tmb$nllfactor,
    randomF = df.tmb$randomF,
    randomM = df.tmb$randomM,
    randomR = df.tmb$randomR,
    recmodel = df.tmb$recmodel,
    nenv = df.tmb$nenv

  )




return(df.sim)


}
