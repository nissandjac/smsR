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
#' fmsy <- getFmsy.OM(OM = om_run, df = om_params, nyears = 100,
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

  nage <- df$nage
  nseason <- df$nseason
  tEnd <- df$nyear # last historical year index

  # Maturity, mortality, weight, selectivity: hold the last historical year's
  # (age x season) pattern constant through the nyears grid-search horizon.
  # Plain array(x, dim=...) recycling does not preserve season alignment once
  # nseason > 1 - it interleaves seasons across years - so .repeat_forecast_year()
  # (defined in forecast_op.R) is used instead.
  mat <- df$mat
  if (is.null(dim(mat)) == FALSE) mat <- mat[, tEnd, , , drop = TRUE]
  M <- df$M
  if (is.null(dim(M)) == FALSE) M <- M[, tEnd, , , drop = TRUE]
  weca <- df$weca
  if (is.null(dim(weca)) == FALSE) weca <- weca[, tEnd, , , drop = TRUE]
  west <- df$west
  if (is.null(dim(west)) == FALSE) west <- west[, tEnd, , , drop = TRUE]

  if (length(dim(mat)) != 4) mat <- .repeat_forecast_year(mat, nyears, nage, nseason)
  if (length(dim(M)) != 4) M <- .repeat_forecast_year(M, nyears, nage, nseason)
  if (length(dim(weca)) != 4) weca <- .repeat_forecast_year(weca, nyears, nage, nseason)
  if (length(dim(west)) != 4) west <- .repeat_forecast_year(west, nyears, nage, nseason)

  propM <- array(.repeat_forecast_year(df$propM[, tEnd, ], nyears, nage, nseason), dim = c(nage, nyears, nseason))
  propF <- array(.repeat_forecast_year(df$propF[, tEnd, ], nyears, nage, nseason), dim = c(nage, nyears, nseason))

  # Selectivity may already be a full time-indexed array or a plain age vector.
  Fsel <- df$Fsel
  if (is.null(dim(Fsel)) == FALSE) Fsel <- Fsel[, tEnd, , , drop = TRUE]
  if (length(dim(Fsel)) != 4) Fsel <- .repeat_forecast_year(Fsel, nyears, nage, nseason)
  Fsel <- array(Fsel, dim = c(nage, nyears, nseason)) # drop the singleton space dim

  Ninit <- df$Ninit

  # Mirror the recruitment models run.agebased.sms.op() (via simulate_om_dynamics())
  # actually supports, deriving each one's required parameters from df.
  recruitment_types <- c("Ricker", "hockey", "Rmean", "estimated", "BH_steep", "BH_env")
  if (!(recruitment %in% recruitment_types)) {
    stop("recruitment must be one of: ", paste(recruitment_types, collapse = ", "))
  }

  alphaSR <- NA
  betaSR <- NA
  h_val <- NA
  Rin <- NA
  R0 <- df$R0

  if (recruitment == "hockey") {
    betaSR <- df$betaSR
    alphaSR <- df$R0 / df$betaSR
  } else if (recruitment == "Ricker") {
    alphaSR <- df$alpha
    betaSR <- df$betaSR
  } else if (recruitment == "estimated") {
    # There is no future df$Rin to draw on for an equilibrium run; hold
    # recruitment at the geometric mean of the historical series for the
    # whole grid-search horizon (same idea used in forecast_op()).
    Rin <- rep(exp(mean(log(df$Rin))), nyears)
  } else if (recruitment %in% c("BH_steep", "BH_env")) {
    h_val <- df$h
  }

  if (!is.null(df$env)) {
    beta_env <- df$beta_env

    # For environmental input, repeat the last historical year
    env <- matrix(df$env[tEnd, ], ncol = ncol(df$env), nrow = nyears, byrow = TRUE)
  } else {
    env <- NULL
    beta_env <- NULL
  }

  # Survey catchability is already available on df, not year-indexed.
  Q <- df$Q
  Q[is.na(Q)] <- 0



  df.fmsy <- list(
    years = years,
    nseason = nseason,
    age = df$age,
    nage = nage,
    F0 = matrix(0, length(years), nseason),
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
    rseason = df$rseason,
    Fmodel = "sim",
    Ninit = Ninit,
    Rin = Rin,
    move = FALSE,
    R0 = R0,
    h = h_val,
    logSDR = log(0),
    logSDcatch = log(0),
    beta_env = beta_env,
    env = env,
    alpha = alphaSR,
    b = rep(0, length(years))
  )


  tmp0 <- run.agebased.sms.op(df.fmsy)

  Fmsyin <- seq(0.01 / max(Fsel), 2 / max(Fsel), length.out = 50) # scale to selectivity

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
