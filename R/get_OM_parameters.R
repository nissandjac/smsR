#' Prepare parameters for an operating model using a fitted sms model
#' @param df.tmb list of smsR input data
#' @param sas fitted smsR stock assessment
#' @param surveySD SD on survey
#' @param recruitment Type of recruitment function shape used by the returned OM
#'   (consumed by \code{\link{forecast_op}} / \code{\link{run.agebased.sms.op}}).
#'   Defaults to \code{NULL}, which resolves to \code{"estimated"} (recruitment
#'   held at the geometric mean of the historical series, no SSB feedback) -
#'   the usual choice for short-term forecasts even when a reference point like
#'   Fmsy was derived under a stock-recruit assumption. Pass \code{"hockey"}
#'   explicitly to instead forecast under a hockey-stick stock-recruit
#'   relationship; when \code{df.tmb$recmodel == 2} (recruitment estimated
#'   freely per year in the assessment, no analytic SR relationship), that
#'   hockey-stick is fit post-hoc to the assessment's own (SSB, R) history via
#'   \code{method}, identically to \code{\link{getFmsy}}, so a forecast run at
#'   \code{Fnew = getFmsy(..., recruitment = "hockey")$Fmsy} reproduces the
#'   same equilibrium dynamics Fmsy was computed under.
#' @param method Only used when \code{recruitment == "hockey"} and
#'   \code{df.tmb$recmodel == 2}. \code{"nls"} (default) or \code{"segmented"} -
#'   see \code{\link{getFmsy}}.
#' @param nspace number of spatial cells
#' @param moveinit Initial distribution in spatial cells (must add to 1)
#' @param movemax maximum movement rate per cell
#' @param rec.space relative recruitment in cells
#' @param moveslope slope of the movement function shape
#' @param movefifty age at 50 percent movement rate out of a cell
#'
#' @return a list of parameters for the OM
#' @export
#'
#' @importFrom abind abind
#'
get_OM_parameters <- function(df.tmb,
                              sas = NULL,
                              surveySD = 0.4,
                              recruitment = NULL,
                              method = "nls",
                              nspace = 1,
                              moveinit = 1,
                              movemax = 0.3,
                              rec.space = 1,
                              moveslope = .7,
                              movefifty = 1) {
  # Do the movement parameters


  if (length(moveinit) == 1) {
    moveinit <- rep(1 / nspace, nspace)
  }
  # Maturity

  # Assign movement out of area
  if (length(movemax) == 1) {
    movemax <- rep(movemax, df.tmb$nseason)
  }

  # if(length(movemax) != nspace){
  #   stop('insert movement rates for each area')
  # }


  movemat <- array(0, dim = c(df.tmb$nage, df.tmb$nyear, nspace, df.tmb$nseason)) # Chances of moving in to the other grid cell
  age <- df.tmb$age


  if (nspace == 1) {
    move <- FALSE
  } else {
    move <- TRUE
  }

  if (move == TRUE) {
    for (j in 1:nspace) {
      for (i in 1:df.tmb$nseason) {
        movemat[, , j, i] <- movemax[j] / (1 + exp(-moveslope * (age - movefifty)))
      }
    }
    movemat[1, , , ] <- 0 # Recruits and 1 year olds don't move
  }


  if (is.null(sas) == FALSE) {
    parms.true <- getEstimatedParms(sas)


    if (is.null(recruitment)) {
      rec <- exp(parms.true$value[parms.true$parameter == "logRin"])
    }


    if(df.tmb$recmodel == 2){
      rec <- getR(df.tmb, sas)$R
    }

  }

  # Recruitment shape for the returned OM. Defaults to "estimated" (geometric
  # mean of historical recruitment, no SSB feedback) regardless of recmodel -
  # short-term forecasts commonly use this even when Fmsy was derived under a
  # stock-recruit assumption. Previously this default was hardcoded into the
  # returned list itself, so passing recruitment = "hockey" silently had no
  # effect; now the argument is actually honored. Pass recruitment = "hockey"
  # explicitly to get a forecast that matches getFmsy(..., recruitment = "hockey").
  if (is.null(recruitment)) {
    recruitment <- "estimated"
  }

  alphaSR <- NA
  betaSR  <- df.tmb$betaSR
  R0      <- df.tmb$betaSR * exp(parms.true$value[parms.true$parameter == "logalpha"])

  if (recruitment == "hockey") {
    if (df.tmb$recmodel == 1) {
      # recmodel == 1 fits an analytic hockey-stick inside the assessment
      # itself, so logalpha is already the right estimate.
      alphaSR <- exp(parms.true$value[parms.true$parameter == "logalpha"])
      betaSR  <- df.tmb$betaSR
    } else {
      # recmodel == 2 has no analytic SR relationship - fit one post-hoc to
      # the assessment's (SSB, R) history, exactly as getFmsy() does, so the
      # two stay consistent.
      hs_fit  <- .fit_hockey_recruitment(df.tmb, sas, method = method)
      alphaSR <- hs_fit$alphaSR
      betaSR  <- hs_fit$betaSR
    }
    R0 <- alphaSR * betaSR
  }



  # Turn life history parameters into spatial objects

  if(is.null(sas)){
    stop('add smsR object to calculate fishing mortality and selectivty')
  }
  F0 <- getF(df.tmb, sas)
  Fsel <- getSel(df.tmb, sas)

  # Into matrix
  F0_flat <- array(F0$F0, dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))
  mat_flat <- array(as.numeric(df.tmb$Mat[, 1:df.tmb$nyears, ]), dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))
  weca_flat <- array(as.numeric(df.tmb$weca[, 1:df.tmb$nyears, ]), dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))
  west_flat <- array(as.numeric(df.tmb$west[, 1:df.tmb$nyears, ]), dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))
  M_flat <- array(as.numeric(df.tmb$M[, 1:df.tmb$nyears, ]), dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))
  Fsel_flat <- array(Fsel$Fsel, dim = c(df.tmb$nage, df.tmb$nyears, 1, df.tmb$nseason))

  # Abind to two spatial objects
  # This assumes the same M, weca, F, and mat in the number of areas
  if(nspace > 1){
    for (i in 1:(nspace - 1)) {
      if (i == 1) {
        F0 <- F0_flat
        mat <- mat_flat
        weca <- weca_flat
        west <- west_flat
        M <- M_flat
        Fsel <- Fsel_flat
      }else{

        F0 <- abind::abind(F0, F0_flat, along = 3)
        mat <- abind::abind(mat, mat_flat, along = 3)
        west <- abind::abind(west, west_flat, along = 3)
        weca <- abind::abind(weca, weca_flat, along = 3)
        M <- abind::abind(M, M_flat, along = 3)
        Fsel <- abind::abind(Fsel, Fsel_flat, along = 3)
      }
    }
  }else{
    F0 <- F0_flat
    mat <- mat_flat
    weca <- weca_flat
    west <- west_flat
    M <- M_flat
    Fsel <- Fsel_flat
  }

  Q <- getCatchability(df.tmb, sas)
  Q <- array(Q$Q, dim = c(df.tmb$nage, df.tmb$nsurvey), )
  Q[is.na(Q)] <- 0






  df.OM <- list(
    years = df.tmb$years,
    nseason = df.tmb$nseason,
    nspace = nspace,
    nyear = df.tmb$nyears,
    movemat = movemat,
    age = df.tmb$age,
    nage = length(df.tmb$age),
    F0 = F0,
    Catchobs = df.tmb$Catchobs,
    Surveyobs = df.tmb$Surveyobs, # Historical survey and catch might be observed values
    M = M,
    mat = mat,
    weca = weca,
    west = west,
    Fsel = Fsel,
    propF = df.tmb$propF,
    propM = df.tmb$propM,
    Fbarage = df.tmb$Fbarage,
    betaSR = betaSR,
    nsurvey = df.tmb$nsurvey,
    surveyStart = df.tmb$surveyStart,
    surveyEnd = df.tmb$surveyEnd,
    surveySD = surveySD,
    surveySeason = df.tmb$surveySeason,
    Q = Q,
    recruitment = recruitment,
    rec.space = rec.space,
    rec.model = df.tmb$recmodel,
    rseason = df.tmb$recseason,
    Fmodel = "est",
    Ninit = c(
      0,
      exp(parms.true$value[parms.true$parameter == "logNinit"])
    ),
    Rin = rec,
    move = move,
    R0 = R0,
    alpha = alphaSR,
    logSDR = parms.true$value[parms.true$parameter == "logSDrec"],
    logSDcatch = log(0),
    b = rep(0, df.tmb$nyears),
    last_year = max(df.tmb$years)
  )


  return(df.OM)
}
