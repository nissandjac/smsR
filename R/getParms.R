#' Title
#'
#' @param df.tmb input values for TMB file
#'
#' @return
#' @export
#'
#' @examples
getParms <- function(df.tmb){


  # Try initial parameters in weird spots
  SDCatch <- rep(.6,max(df.tmb$Cidx_CV)+1)
  SDsurvey <- rep(.1, length(unique(df.tmb$Qidx_CV[df.tmb$Qidx_CV > -1])))
  logQ <- sum(df.tmb$Qlastage-df.tmb$Qminage)+length(df.tmb$Qidx)

  parms <- list(logRin = rep(18, df.tmb$nyears),
                  logNinit = rep(15, df.tmb$nage-1),
                  Fyear = rep(1, df.tmb$nyears), # Mapped out
                  Fseason = matrix(0.5, nrow = 1, ncol = length(unique(df.tmb$bidx))),
                  logFage = matrix(1, nrow = df.tmb$Fmaxage+1, ncol = length(unique(df.tmb$bidx))),
                  SDsurvey = SDsurvey,
                  SDcatch = as.matrix(SDCatch),
                  logQ = rep(log(1), logQ),#length(df.tmb$surveyCV)
                  pin = 1,
                  logalpha = 2,
                  logbeta = df.tmb$betaSR,
                  logSDrec = log(0.5))





  return(parms)
}



#' Title
#'
#' @param df.tmb list of tmb input parameters
#'
#' @return
#' @export
#'
#' @examples
getMPS <- function(df.tmb, parms, mapExtra = NA){

  # Not estimated parameters
  mps <- list()

  if(sum(df.tmb$power) == 0){
    mps$pin <- factor(parms$pin*NA)
  }


  if(df.tmb$estCV[2] == 2){
    mps$SDcatch <- factor(parms$SDcatch*NA)
  }

  if(df.tmb$useEffort == 1){
    mps$Fyear <- factor(parms$Fyear*NA)
  }

  if(is.numeric(df.tmb$betaSR)){
    mps$logbeta <- factor(parms$logbeta*NA)
  }

  for(i in 1:length(mapExtra)){

   if(is.na(mapExtra[i]) ==0){

       mps[[mapExtra[i]]] <- factor(parms[[mapExtra[i]]]*NA)
     }

   }


return(mps)

}


