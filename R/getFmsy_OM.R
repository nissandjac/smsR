#' Estimate Fmsy from a SMS Operating Model
#'
#' Calculates the fishing mortality rate (\emph{F}) that produces maximum sustainable yield (Fmsy)
#' using an SMS operating model. This function simulates forward under a
#' specified recruitment assumption and identifies the \emph{F} that maximizes yield.
#'
#' @param OM An sms operating model from \code{"run.agebased.sms.op()"}
#' @param df A list containing operating model parameters.
#' @param nyears Integer. Number of years to simulate (default is 50).
#' @param recruitment Character. Recruitment model to use in the simulation. Currently only \code{"hockey"} is supported.
#' @param plotMSY Logical. If \code{TRUE}, a plot of the yield curve and Fmsy estimate is produced.
#' @param stochastic Logical. Should stochastic recruitment deviations be included in the simulation?
#'
#' @return A numeric value representing the fishing mortality rate that results in maximum sustainable yield (Fmsy).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fmsy <- getFmsy.OM(sas = fitted_model, df = om_params, nyears = 100,
#'                    recruitment = "hockey", plotMSY = TRUE, stochastic = TRUE)
#' }
getFmsy.OM <- function(OM,
                       df,
                    nyears = 50,
                    recruitment = "hockey",
                    plotMSY = FALSE,
                    stochastic = FALSE) {
  #  Number
  years <- 1:nyears
  nyears <- length(years)


  # Maturity
  mat <- df$mat # Fecundity
  M <- df$M
  weca <- df$weca
  west <- df$west
  #
  if (length(dim(M)) != 4) {
    M <- array(M, dim = c(df$nage, df$nyear + 1, 1, df$nseason))
  }

  if (length(dim(weca)) != 4) {
    weca <- array(weca, dim = c(df$nage, df$nyear + 1, 1, df$nseason))
  }

  if (length(dim(west)) != 4) {
    west <- array(west, dim = c(df$nage, df$nyear + 1, 1, df$nseason))
  }

  if (length(dim(mat)) != 4) {
    mat <- array(mat, dim = c(df$nage, df$nyear + 1, 1, df$nseason))
  }

  M <- array(M[, df$nyear, ], dim = c(df$nage, nyears, df$nseason))
  mat <- array(Mat[, df$nyear, ], dim = c(df$nage, nyears, df$nseason))
  weca <- array(weca[, df$nyears, ], dim = c(df$nage, nyears, df$nseason))
  west <- array(west[, df$nyears, ], dim = c(df$nage, nyears, df$nseason))
  propM <- array(df$propM[, df$nyears, ], dim = c(df$nage, nyears, df$nseason))
  propF <- array(df$propF[, df$nyears, ], dim = c(df$nage, nyears, df$nseason))

  Fselin <- getSelectivity(df, sas)
  #Fselin <- Fselin[Fselin$block == max(Fselin$block),]
  # I like this scaled to 1 for consistency

  # blocks <- unique(Fselin$blocks)


  # for(i in 1:length(blocks)){
  #
  #   Fsel$Fsel[Fsel$blocks == blocks[i]] <- Fsel$Fsel[Fsel$blocks == blocks[i]]/max(Fsel$Fsel[Fsel$blocks == blocks[i]])
  #
  # }

  for (i in 1:df$nseason) {
    if (i == 1) {
      Fsel <- array(Fselin$selec[Fselin$season == i], dim = c(df$nage, nyears, 1))
    } else {
      tmp <- array(Fselin$selec[Fselin$season == i], dim = c(df$nage, nyears, 1))
      Fsel <- abind::abind(Fsel, tmp, along = 3)
    }
  }

  parms.est <- getEstimatedParms(sas)

  Ninit <- exp(parms.est$value[parms.est$parameter == "logNinit"])
  # SDR <- exp(parms.est$value[parms.est$parameter == 'logSDrec'])

  if (recruitment == "hockey") {
    alphaSR <- exp(parms.est$value[parms.est$parameter == "logalpha"])
    betaSR <- df$betaSR
    R0 <- alphaSR * df$betaSR
  } else {
    alphaSR <- NA
    betaSR <- NA
    R0 <- exp(parms.est$value[parms.est$parameter == "logR0"])
  }

  if(sas$dat$nenv >0){

    beta_env <- sas$reps$par.fixed[names(sas$reps$par.fixed) == 'env']

    # For environmental input, repeat the last year
    env <- matrix(df$env_matrix[df$nyears], ncol = df$nenv, nrow = nyears)

  }else{
    env <- NULL
    beta_env <- NULL
  }

  # Survey catchability
  Q <- getCatchability(df, sas)
  Q <- array(Q$Q, dim = c(df$nage, df$nsurvey), )
  Q[is.na(Q)] <- 0



  df.fmsy <- list(
    years = years,
    nseason = df$nseason,
    age = df$age,
    nage = length(df$age),
    F0 = matrix(0, length(years), df$nseason),
    Fsel = Fsel,
    propM = propM,
    propF = propF,
    M = M,
    mat = mat,
    weca = weca,
    west = west,
    sel = Fsel,
    betaSR = betaSR,
    Fbarage = df$Fbarage,
    nsurvey = df$nsurvey,
    surveySeason = df$surveySeason,
    surveyStart = df$surveyStart,
    surveyEnd = df$surveyEnd,
    surveySD = 0,
    Q = Q,
    recruitment = recruitment,
    rseason = df$recseason,
    Fmodel = "sim",
    Ninit = NULL,
    Rin = NA,
    move = FALSE,
    R0 = R0,
    h = exp(-0.6931472), # Fix this later. need to find it in 'sas'
    logSDR = log(0),
    beta_env = beta_env,
    env = env,
    alpha = alphaSR,
    b = rep(0, length(years))
  )


  tmp0 <- run.agebased.sms.op(df.fmsy)

  Fmsyin <- seq(0.01 / max(Fselin$selec), 2 / max(Fselin$selec), length.out = 50) # scale to selectivity

  for (i in 1:length(Fmsyin)) {
    df.fmsy$F0 <- matrix(Fmsyin[i], length(years), df$nseason)

    tmp <- run.agebased.sms.op(df.fmsy)

    if (i == 1) {
      df.out <- data.frame(F0 = tmp$Fbar[length(tmp$Fbar)], catch = tmp$Catch[length(years)], SSB = tmp$SSB[length(years)])
    } else {
      df.out <- rbind(
        df.out,
        data.frame(F0 = tmp$Fbar[length(tmp$Fbar)], catch = tmp$Catch[length(years)], SSB = tmp$SSB[length(years)])
      )
    }
  }

  if (plotMSY == TRUE) {
    print(ggplot(df.out, aes(x = F0, y = catch)) +
            geom_line() +
            theme_classic() +
            scale_x_continuous("Fbar"))
  }


  Fmsy <- list(
    Fmsy = df.out$F0[which.max(df.out$catch)],
    MSY = max(df.out$catch)
  )


  return(Fmsy)
}
