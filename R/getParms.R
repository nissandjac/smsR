#' Title
#'
#' @param df.tmb input values for TMB file
#' @param parms.true data frame of true parameters
#'
#' @return
#' @export
#'
#' @examples
getParms <- function(df.tmb, parms.true = NULL){


  # Try initial parameters in weird spots
    SDCatch <- rep(.6,max(df.tmb$Cidx_CV)+1)
    SDsurvey <- rep(.1, length(unique(df.tmb$Qidx_CV[df.tmb$Qidx_CV > -1])))
    logQ <- sum(df.tmb$Qlastage-df.tmb$Qminage)+length(df.tmb$Qidx)

    if(is.null(df.tmb$betaSR)){
      betaSR <-50000
    }else{
      betaSR <- df.tmb$betaSR
    }

    parms <- list(logRin = rep(18, df.tmb$nyears),
                  logNinit = rep(15, df.tmb$nage-1),
                  Fyear = rep(1, df.tmb$nyears), # Mapped out
                  Fseason = matrix(1, nrow = df.tmb$nseason-1, ncol = length(unique(df.tmb$bidx))),
                  logFage = matrix(1, nrow = length(df.tmb$Fminage:df.tmb$Fmaxage), ncol = length(unique(df.tmb$bidx))),
                  SDsurvey = SDsurvey,
                  SDcatch = as.matrix(SDCatch),
                  logQ = rep(log(1), logQ),#length(df.tmb$surveyCV)
                  pin = 1,
                  logalpha = 2,
                  logbeta = log(betaSR),
                  logSDrec = log(0.5))



    if(!is.null(parms.true)){

      pnames <- unique(parms.true$parameter)


      for(i in 1:length(pnames)){

        if(is.null(dim(parms[[i]]))){
        parms[names(parms) == pnames[i]][[1]] <- parms.true$value[parms.true$parameter == pnames[i]]

        }else{

        if(length(dim(parms[[i]])) == 2){
          parms[names(parms) == pnames[i]] <- array(parms.true$value[parms.true$parameter == pnames[i]],
                                                    dim = c(dim(parms[[i]])[1],dim(parms[[i]])[2]
                                                    )
                                                    )
        }

        if(length(dim(parms[[i]])) == 3){
        parms[names(parms) == pnames[i]] <- array(parms.true$value[parms.true$parameter == pnames[i]],
                                                  dim = c(dim(parms[[i]])[1],dim(parms[[i]])[2],
                                                          dim(parms[[i]])[3]
                                                  )
        )
        }

        }

      }
    }


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

  if(df.tmb$nseason == 1){
    mps$Fseason <- factor(parms$Fseason*NA)
  }


  for(i in 1:length(mapExtra)){

   if(is.na(mapExtra[i]) ==0){

       mps[[mapExtra[i]]] <- factor(parms[[mapExtra[i]]]*NA)
     }

   }


return(mps)

}


