#' Ready parameters for an operating model using a fitted sms model
#'
#' @param df.tmb
#' @param sas
#' @param move
#' @param surveySD
#' @param recruitment
#'
#' @return
#' @export
#'
#' @examples
get_OM_parameters <- function(df.tmb, sas = NULL,
                              surveySD = 0.4,
                              recruitment = NULL,
                              nspace = 1,
                              moveinit = 1,
                              movemaxinit = 0.3,
                              moveslope = .7,
                              movefifty = 1
                              ){


   # Do the movement parameters


  if(moveinit ==1){
    moveinit <- rep(1/nspace, nspace)
  }
  # Maturity

  # Assign movement out of area
  movemax <- rep(movemaxinit,nseason)

  movemat <- array(0, dim = c(nspace, df.tmb$nage, df.tmb$nseason, df.tmb$nyear)) # Chances of moving in to the other grid cell
  age <- df.tmb$age


  if(nspace == 1){
    move = FALSE
  }else{
    move = TRUE
  }

  if(move == TRUE){
    for(j in 1:nspace){
      for(i in 1:nseason){
        movemat[j,,i,] <- movemax[i]/(1+exp(-moveslope*(age-movefifty)))

      }
    }
    movemat[,1:2,,] <- 0 # Recruits and 1 year olds don't move

  }


    if(is.null(sas) == FALSE){
      parms.true <- getEstimatedParms(sas)
    }

    if(is.null(recruitment)){
      rec <- exp(parms.true$value[parms.true$parameter == 'logRin'])
    }

    # Turn life history parameters into spatial objects
    F0 <- getF(df.tmb, sas)
    # Into matrix
    F0 <- array(F0$F0, dim = c(df.tmb$nage,df.tmb$nyears,1,df.tmb$nseason))
    mat <- df.tmb$Mat
    weca <- df.tmb$weca
    west <- df.tmb$west
    # Abind to two spatial objects
    for(i in 1:nspace){

    }


    df.OM <- list(
    years = df.tmb$years,
    nseason = df.tmb$nseason,
    age = df.tmb$age,
    nage = length(df.tmb$age),
    F0 = df.tmb$f0,
    sel = mtrx$Fsel,
    M = mtrx$M,
    mat = mtrx$mat,
    weca = mtrx$weca,
    west = mtrx$west,
    propF = df.tmb$propF,
    propM = df.tmb$propM,
    Fbarage = df.tmb$Fbarage,
    betaSR = df.tmb$betaSR,
    nsurvey = df.tmb$nsurvey,
    surveyStart = df.tmb$surveyStart,
    surveyEnd = df.tmb$surveyEnd,
    surveySD = surveySD,
    surveySeason = df.tmb$surveySeason,
    Q = Q,
    recruitment = 'estimated',
    rseason = df.tmb$recseason,
    Fmodel = 'est',
    Ninit = c(0,
              exp(parms.true$value[parms.true$parameter == 'logNinit'])),
    Rin = rec,
    move = move,
    R0 =  df.tmb$betaSR*exp(parms.true$value[parms.true$parameter == 'logalpha']),
    SDR = exp(parms.true$value[parms.true$parameter == 'logSDrec']),
    b = rep(0, df.tmb$nyears)
  )






return(df.OM)

}
