#' Run a full age based model as operating model for SMS
#'
#' @param df list of parameters and life history values
#'
#' @return
#' returns a list of derived and calculated values
#' @export
#'
#' @examples
#' run.agebased.true.catch(df) # runs the model
#'
run.agebased.sms.op <- function(df,
                                seed = NULL,
                                OM = NULL) {

  if(is.null(seed) == 1){
    set.seed(round(runif(1, min = 1, max = 1e6),digits = 0))
  }else{
    set.seed(seed)
  }

  nseason <- df$nseason

  df$tEnd <- length(df$years) * nseason
  nyear <- df$tEnd / df$nseason #
  year <- df$years

  # Maturity
  mat <- df$mat # Fecundity
  mat[is.na(mat)] <- 0


  # Mortality
  M <- df$M
  M[is.na(M)] <- 0 # Remove NA's

  # Weight
  weca <- df$weca
  weca[is.na(weca)] <- 0

  west <- df$west
  west[is.na(west)] <- 0

  # Fishing mortality
  F0 <- df$F0
  F0[is.na(F0)] <- 0
  # Age
  nage <- df$nage
  age <- df$age
  maxage <- max(df$age)

  R0 <- df$R0

  # Calculate N0 based on R0
  mage <- max(df$age) # Max age
  agetmp <- 0:(mage * 3)
  nagetmp <- length(agetmp)

  if (length(dim(F0)) != 4) {
    F0 <- array(F0, dim = c(nage, nyear, 1, nseason))
  }

  selin <- df$Fsel

  if (length(dim(selin)) != 4) {
    selin <- array(selin, dim = c(nage, nyear, 1, nseason))
  }


  if (is.null(df$nspace)) {
    df$nspace <- 1
  }


  if(length(dim(M)) != 4){
    M <- array(M, dim = c(df$nage, ncol(M), 1, df$nseason))
  }


  M0 <- rep(0, mage * 3 + 1)
  M0[1:nage] <- M[, 1, 1, 1] * df$nseason
  M0[nage:length(M0)] <- M[nage, 1, 1, 1] * df$nseason

  if(length(dim(west)) != 4){ # add a one to the spatial dimension
    west <- array(west, dim = c(df$nage, ncol(west), 1, df$nseason))
  }

  if(length(dim(mat)) != 4){ # Same as above
    mat <- array(mat, dim = c(df$nage, ncol(mat), 1, df$nseason))
  }

  if(length(dim(weca)) != 4){ # Same as above
    weca <- array(weca, dim = c(df$nage, ncol(weca), 1, df$nseason))
  }


  if(df$recruitment == 'BH_steep'){
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

  # Initialize for saving
  year_1 <- c(year, max(year) + 1)

  # Distribute over space
  if (is.null(df$Ninit)) {
    Ninit <- N0
  } else {
    Ninit <- df$Ninit
  }

  if(is.null(df$rec.space)){
      df$rec.space <- rep(1/df$nspace, df$nspace)
  }

  if (is.null(df$movemat)) {
    movemat <- array(0, dim = c(df$nage, nyear, 1, df$nseason))
  } else {
    movemat <- df$movemat
  }

  # Fishing mortality is extracted from df: 'est' uses the estimated F-at-age
  # array directly, 'sim' rebuilds it from a simulated F trend and selectivity.
  if (df$Fmodel == "est") {
    Fseason_fun <- function(yr, season, space, sel_season) {
      F0[, yr, space, season]
    }
  } else if (df$Fmodel == "sim") {
    Fseason_fun <- function(yr, season, space, sel_season) {
      df$F0[yr, season] * sel_season
    }
  } else {
    stop("Unknown df$Fmodel: must be 'est' or 'sim'")
  }

  # Historical run: recruitment for 'estimated' is simply the model-estimated
  # value for that year.
  Rin_fun <- function(yr, space) {
    df$Rin[yr] * df$rec.space[space]
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
    env = df$env,
    stochastic = TRUE
  )

  # Add error to the survey for the assessment model
  eps <- rnorm(length(sim$survey), mean = 0, sd = df$surveySD)
  eps <- array(eps, dim = dim(sim$survey))
  # True survey
  survey.true <- sim$survey
  survey <- sim$survey * exp(eps - 0.5 * df$surveySD^2)
  survey[survey < 0] <- -1

  # Add Error to catches
  eps_c <- rnorm(sim$Catch.save.age, mean = 0, sd = exp(df$logSDcatch))
  eps_c <- array(eps_c, dim = dim(sim$Catch.save.age))
  Catchobs <- (sim$CatchN.save.age * exp(eps_c - 0.5 * exp(df$logSDcatch)^2))
  Catchobs <- apply(Catchobs, MARGIN = c(1,2,4), sum)

  df.out <- list(
    N.save = sim$N.save,
    SSB = sim$SSB,
    SSB0 = SSB_0,
    N.save.age = sim$N.save.age,
    R.save = sim$R.save,
    R.err.save = sim$R.err.save,
    V.save = sim$V.save,
    E.save = sim$E.save,
    SSB.all = sim$SSB.all,
    Catch.save.age = sim$Catch.save.age,
    CatchN.save.age = sim$CatchN.save.age,
    Catchobs = Catchobs,
    Catch = sim$Catch,
    Catch.age = sim$Catch.age,
    survey = survey,
    survey.true = survey.true,
    Fbar = sim$Fbar,
    M = sim$M,
    age_comps_OM = sim$age_comps_OM,
    age_catch = sim$age_catch,
    Z = sim$Z,
    Fseason = sim$Fseason
  )

  return(structure(df.out, class = 'sms_om'))
}
