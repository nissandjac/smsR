#' Title
#'
#' @param OM
#' @param df
#' @param recmodel
#' @param MCV
#' @param M_max
#'
#' @return
#' @export
#'
#' @examples
OM_to_smsR <- function(OM,
                      df,
                      recmodel = 1,
                      MCV = NULL,
                      M_max = max(df$age)){

  #
  nyear <- df$nyear
  nage <- df$nage
  nseason <- df$nseason

  Surveyobs <- OM$survey
  Catchobs <- array(OM$CatchN.save.age[,,,1], dim = c(df$nage, nyear,  df$nseason)) # smsR doesn't handle the spatial dimension
  Catchobs_err <- Catchobs
  # Add uncertainty to catch
  for(i in 1:nyear){
  Catchobs_err[,i,] <- Catchobs[,i,] * exp(rnorm(nage, mean = 0, sd = exp(df$logSDcatch)))
  }


  # Add uncertainty to survey
  Surveyobs[Surveyobs == -1] <- NA
  for(i in 1:nyear){
    for(j in 1:nage){
      for(k in 1:df$nsurvey){
        Surveyobs[j,i,k] <- Surveyobs[j,i,k] * exp(rnorm(1, mean = 0, sd = df$surveySD))
     }
    }
  }
  Surveyobs[is.na(Surveyobs)] <- -99



  # Group index and effort #
  weca <- df$weca
  west <- df$west
  M <- df$M
  mat <- df$mat

  if (length(dim(M)) != 4) {

    if(is.null(dim(M))){
    M <- c(M, M[nyear])  # replace with your actual vector

    # Step 2: Repeat across 11 rows
    M <- matrix(rep(M, each = nage), nrow = nage, ncol = nyear+1)

    # Step 3: Convert to 3D array (11, 31, 1)
    M <- array(M, dim = c(nage, nyear + 1 , nseason))
    }else{
    M <- array(rep(M, nyear+ 1), dim = c(nage, nyear + 1, nseason))
    }
  }

  if (length(dim(weca)) != 4) {
    weca <- array(weca, dim = c(nage, nyear + 1, nseason))
  }

  if (length(dim(west)) != 4) {
    west <- array(west, dim = c(nage, nyear + 1, nseason))
  }

  if (length(dim(mat)) != 4) {
    mat <- array(mat, dim = c(nage, nyear + 1, nseason))
  }

  mtrx <- list(west = west, weca = weca,
               mat = mat,
               M = M)

  # find the minimum catch age


  # Calculate lists for survey and catch CV

  surveySD <- list()

  for(i in 1:df$nsurvey){

    surveySD[[i]] <- df$Qminage[i]

  }

  catchSD <- list(df$Fminage)


  # Get the actual Fmaxage
  Fmaxage <- min(df$age[df$Fsel > .95])
  Qlastage <- min(df$age[df$Fsel > .95])

  betaSR <- min(OM$SSB[OM$R.save > median(OM$R.save)]) # lowest SSB with R > median

  #
  if(is.null(MCV)){
    randomM <- 0
  }else{
    randomM <- 1
  }


  dat <- get_TMB_parameters(mtrx = mtrx ,
                            Surveyobs = Surveyobs,
                            Catchobs = Catchobs_err,
                            years = df$years,
                            nseason = df$nseason,
                            nsurvey = df$nsurvey,
                            ages = df$age,
                            Fbarage = df$Fbarage,
                            recseason = df$rseason,
                            Fminage = df$Fminage,
                            Fmaxage = Fmaxage,
                            Qminage = df$Qminage,
                            Qmaxage = df$Qmaxage,
                            Qlastage = df$Qmaxage,
                            minSDsurvey = df$surveySD*.2,
                            minSDcatch = exp(df$logSDcatch)*.2,
                            surveyStart = df$surveyStart,
                            surveyEnd = df$surveyEnd,
                            estSD = c(0,0,0),
                            MCV = MCV,
                            randomM = randomM,
                            M_max = M_max,
                            recmodel = recmodel,
                            betaSR = betaSR,
                            surveySD = surveySD,
                            catchSD  = catchSD
                            )




  return(dat)

}
