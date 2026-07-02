#' Forecast an smsR Operating Model (OM)
#'
#' Simulates future population dynamics using the operating model for a given number of years.
#' Commonly used in Management Strategy Evaluation (MSE) frameworks.
#'
#' @param df List containing historical OM parameters (e.g., mortality, selectivity, maturity).
#' @param nyear Integer. Number of years to forecast.
#' @param Fnew Numeric or vector. Fishing mortality to apply in the forecast years.
#' @param env_new Optional environmental input(s) used in forecasting (e.g., temperature, productivity).
#'   Row \code{i} is the covariate for forecast year \code{i}.
#' @param Ninit Numeric vector. Numbers-at-age going into the first forecast year, i.e. the
#'   state \emph{after} the last historical year has been survived - the
#'   \code{nyears + 1} slot of \code{run.agebased.sms.op()}'s \code{N.save.age} output
#'   (\code{OM$N.save.age[, dim(OM$N.save.age)[2], , 1]}), not the \code{nyears} slot.
#' @param stochastic Logical. Should recruitment deviations be stochastic? Defaults to \code{FALSE}.
#' @param seed Optional integer. Random seed for reproducibility if \code{stochastic = TRUE}.
#' @param Rmean_years Only used when \code{df$recruitment == "estimated"} (i.e. the
#'   historical run used a free recruitment estimate per year, not a stock-recruit
#'   curve). There is no future \code{df$Rin} to draw on, so forecast recruitment is
#'   instead the geometric mean of the most recent \code{Rmean_years} historical
#'   \code{df$Rin} values (with the usual bias-corrected lognormal deviation applied).
#'   Defaults to \code{NULL}, which averages over the full historical \code{df$Rin} series.
#'   Ignored for stock-recruit-based recruitment models (\code{"Ricker"}, \code{"hockey"},
#'   \code{"BH_steep"}, \code{"BH_env"}) - those already respond to forecast SSB. To
#'   forecast from a stock-recruit relationship instead of a geometric mean, set
#'   \code{df$recruitment} to one of those types before calling \code{forecast_op()}.
#'
#' @return A list containing forecasted operating model outputs, such as numbers-at-age, SSB, recruitment, and catch.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Forecast 10 years with fixed F
#' out <- forecast_om(df = om_data, nyear = 10, Fnew = 0.2, Ninit = N_last_year)
#'
#' # Forecast with stochastic recruitment and environmental input
#' out <- forecast_om(df = om_data, nyear = 20, Fnew = rep(0.2, 20),
#'                    env_new = env_proj, Ninit = N_last_year,
#'                    stochastic = TRUE, seed = 123)
#' }
#'
#' @keywords internal
.repeat_forecast_year <- function(slice, nyear_target, nage, nseason) {
  # Hold a single year's (age x season) pattern constant across the whole
  # forecast horizon. Plain array(x, dim=...) recycling does not preserve
  # season alignment once nseason > 1 - it interleaves seasons across years
  # instead of repeating the same season pattern every year - so this fills
  # season-by-season instead.
  slice <- matrix(slice, nrow = nage, ncol = nseason)
  out <- array(0, dim = c(nage, nyear_target, 1, nseason))
  for (s in 1:nseason) {
    out[, , 1, s] <- slice[, s]
  }
  out
}
forecast_op <- function(df,
                        Ninit,
                        nyear = 2,
                        Fnew = matrix(.2, nrow = df$nseason, ncol = nyear),
                        env_new = NULL,
                        stochastic = FALSE,
                        seed = NULL,
                        Rmean_years = NULL){

  if(is.null(seed) == 1){
    set.seed(round(runif(1, min = 1, max = 1e6),digits = 0))
  }else{
    set.seed(seed)
  }

  if(is.matrix(Fnew) == FALSE){
    Fnew <- matrix(Fnew)
  }

  if(is.null(env_new) == 0){
    if(is.matrix(env_new) == FALSE){
    env_new <- matrix(env_new)
    }
  }


  nseason <- df$nseason
  tEnd <- df$nyear # Last historical year index (M/mat/weca/west are indexed by year, not by season*year)
  # Only 1
  # Make all the matrices two dimensional to save space
  mat <- df$mat # Fecundity
  if(is.null(dim(mat)) == FALSE){
    mat <- mat[,tEnd,,, drop = TRUE]
  }
  # Mortality
  M <- df$M
  if(is.null(dim(M)) == FALSE){
  M <- M[,tEnd,,, drop = TRUE]
  }
  # Weight
  weca <- df$weca[,tEnd,,]
  west <- df$west[,tEnd,,]

  # Age
  nage <- df$nage
  age <- df$age
  maxage <- max(df$age)

  R0 <- df$R0

  # Calculate N0 based on R0
  mage <- max(df$age) # Max age
  agetmp <- 0:(mage * 3)
  nagetmp <- length(agetmp)

  if (length(dim(M)) != 4) {
    M <- .repeat_forecast_year(M, nyear + 1, nage, nseason)
  }

  if (length(dim(weca)) != 4) {
    weca <- .repeat_forecast_year(weca, nyear + 1, nage, nseason)
  }

  if (length(dim(west)) != 4) {
    west <- .repeat_forecast_year(west, nyear + 1, nage, nseason)
  }

  if (length(dim(mat)) != 4) {
    mat <- .repeat_forecast_year(mat, nyear + 1, nage, nseason)
  }

  # Selectivity: df$Fsel may already be a full 4D time series (e.g. from
  # get_OM_parameters()) or a plain age vector with no time dimension (e.g.
  # from sim_OM_parameters()). Either way, hold the last historical year's
  # pattern constant through the forecast.
  selin <- df$Fsel
  if (is.null(dim(selin)) == FALSE) {
    selin <- selin[, tEnd, , , drop = TRUE]
  }

  if (length(dim(selin)) != 4) {
    selin <- .repeat_forecast_year(selin, nyear, nage, nseason)
  }


  if (is.null(df$nspace)) {
    df$nspace <- 1
  }

  if(df$recruitment == 'BH_steep'){
    M0 <- rep(0, mage * 3 + 1)
    M0[1:nage] <- M[, 1, 1, 1] * df$nseason
    M0[nage:length(M0)] <- M[nage, 1, 1, 1] * df$nseason

    N0tmp <- rep(NA, nagetmp)

    N0tmp[1:(nagetmp - 1)] <- R0 * exp(-agetmp[1:(nagetmp - 1)] * M0[1:(nagetmp - 1)])
    N0tmp[nagetmp] <- R0 * exp(-M0[nagetmp] * agetmp[nagetmp]) / (1 - exp(-M0[nagetmp]))

    N0 <- matrix(NA, nage)
    N0[1:(nage - 1)] <- N0tmp[1:(nage - 1)]
    N0[nage] <- sum(N0tmp[nage:nagetmp])

    SSB_0 <- sum(N0 * west[, 1, 1, 1] * mat[, 1, 1, 1])
  }else{
    SSB_0 <- NA
  }

  # Recruitment per country
  nspace <- df$nspace

  # Initialize for saving. Forecast year 1 is the year after the last
  # historical year, not the last historical year itself.
  year <- (df$years[tEnd] + 1):(df$years[tEnd] + nyear)
  year_1 <- c(year, max(year) + 1)

  if(is.null(df$rec.space)){
    df$rec.space <- rep(1/df$nspace, df$nspace)
  }

  # df$propM, df$propF and df$b are read directly by the shared engine and
  # are only sized for the historical period (length/dim = historical nyear).
  # Hold the last historical year's values constant through the forecast so
  # indexing by the forecast's own yr (which can exceed the historical
  # length) doesn't go out of bounds.
  df$propM <- array(.repeat_forecast_year(df$propM[, tEnd, ], nyear, nage, nseason), dim = c(nage, nyear, nseason))
  df$propF <- array(.repeat_forecast_year(df$propF[, tEnd, ], nyear, nage, nseason), dim = c(nage, nyear, nseason))
  df$b <- rep(df$b[tEnd], nyear)

  if (is.null(df$movemat)) {
    movemat <- array(0, dim = c(df$nage, nyear, 1, df$nseason))
  } else {
    # df$movemat is only sized for the historical period; hold the last
    # historical year's movement pattern constant through the forecast
    # (indexing it directly would go out of bounds once nyear exceeds the
    # historical year count).
    nspace_mm <- dim(df$movemat)[3]
    movemat_slice <- df$movemat[, tEnd, , , drop = FALSE]
    movemat <- array(0, dim = c(df$nage, nyear, nspace_mm, df$nseason))
    for (yr in 1:nyear) {
      movemat[, yr, , ] <- movemat_slice[, 1, , ]
    }
  }

  # Fishing mortality in a forecast is whatever the user supplies via Fnew,
  # applied through the (fixed, last-historical-year) selectivity-at-age.
  Fseason_fun <- function(yr, season, space, sel_season) {
    Fnew[1, season] * sel_season
  }

  # When recruitment is 'estimated' there is no future df$Rin to draw on, so
  # forecast recruitment falls back to the geometric mean of recent historical
  # Rin, with the same bias-corrected lognormal deviation used elsewhere. This
  # is unused (and never called) for stock-recruit-based recruitment models.
  if (df$recruitment == "estimated") {
    Rin_hist <- if (is.null(Rmean_years)) df$Rin else utils::tail(df$Rin, Rmean_years)
    R_bar <- exp(mean(log(Rin_hist)))
    SDR_est <- exp(df$logSDR)

    Rin_fun <- function(yr, space) {
      err <- if (stochastic) rnorm(1, mean = 0, sd = SDR_est) else 0
      R_bar * exp(-0.5 * SDR_est^2 + err) * df$rec.space[space]
    }
  } else {
    Rin_fun <- NULL
  }

  sim <- simulate_om_dynamics(
    df = df,
    nyear = nyear,
    year = year,
    year_1 = year_1,
    nage = nage,
    age = age,
    maxage = maxage,
    nseason = nseason,
    nspace = nspace,
    M = M,
    weca = weca,
    west = west,
    mat = mat,
    selin = selin,
    movemat = movemat,
    Ninit = Ninit,
    R0 = R0,
    SSB_0 = SSB_0,
    Fseason_fun = Fseason_fun,
    Rin_fun = Rin_fun,
    env = env_new, # env_new[yr, ] is the covariate for forecast year yr
    stochastic = stochastic
  )

  df.out <- list(
    N.save = sim$N.save,
    SSB = sim$SSB,
    N.save.age = sim$N.save.age,
    R.save = sim$R.save,
    V.save = sim$V.save,
    E.save = sim$E.save,
    SSB.all = sim$SSB.all,
    Catch.save.age = sim$Catch.save.age,
    CatchN.save.age = sim$CatchN.save.age,
    Catch = sim$Catch,
    Catch.age = sim$Catch.age,
    survey = sim$survey,
    survey.true = sim$survey.true,
    Fbar = sim$Fbar,
    age_comps_OM = sim$age_comps_OM,
    age_catch = sim$age_catch,
    Z = sim$Z,
    Fseason = sim$Fseason
  )

  return(df.out)
}
