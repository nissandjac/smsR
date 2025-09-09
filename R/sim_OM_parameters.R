#' Prepare a parameter list for a simulated `smsR` stock (Operating Model)
#'
#' Build a self-contained list of operating-model (OM) inputs for simulating an
#' age-structured stock with optional seasonality, surveys, time-varying fishing
#' mortality, autocorrelated recruitment, and several natural-mortality dynamics.
#'
#' @param nseason Integer. Number of seasons within a year (e.g., 1 = annual).
#' @param nspace Integer. Number of spatial areas (used to split seasonal F if > 1).
#' @param nyear Integer. Number of years to simulate.
#' @param Linf Numeric or vector. Asymptotic length for von Bertalanffy growth.
#'   If numeric, weights/lengths are generated internally; if not, provide
#'   external weight inputs via \code{weca}/\code{west}.
#' @param weca Numeric array/matrix or \code{NULL}. Weight-at-age used for catch.
#'   If \code{NULL} and \code{Linf} is numeric, generated internally.
#' @param west Numeric array/matrix or \code{NULL}. Weight-at-age used for SSB.
#'   If \code{NULL} and \code{Linf} is numeric, generated internally.
#' @param K Numeric. von Bertalanffy growth parameter.
#' @param Fpast Numeric (scalar). Historical constant fishing mortality used only
#'   to shape the initial numbers-at-age distribution (\code{Ninit}). Set 0 to skip.
#' @param t0 Numeric. von Bertalanffy theoretical age at zero length.
#' @param gamma Numeric. Logistic maturity curve slope.
#' @param tau Numeric. Logistic maturity curve inflection age.
#' @param gamma_sel Numeric. Logistic fishing selectivity slope.
#' @param tau_sel Numeric. Logistic fishing selectivity inflection age.
#' @param maxage Integer. Maximum modeled age (age classes are \code{0:maxage}).
#' @param a Numeric. Length–weight coefficient (when generating weights from VBGF).
#' @param b_size Numeric. Length–weight exponent (when generating weights).
#' @param M Numeric. Baseline natural mortality used as the starting level or
#'   constant value depending on \code{mortality}.
#' @param b Numeric or length-\code{nyear} vector. Power for production/biomass
#'   scaling used by the simulator (passed through as \code{b} in the OM list).
#' @param R0 Numeric. Unfished recruitment level (scales recruitment).
#' @param theta Numeric. Observation-model parameter (passed through).
#' @param h Numeric. Beverton–Holt steepness; required when
#'   \code{recruitment \%in\% c("BH_steep","BH_env")}.
#' @param SDcatch Numeric (sd on log scale). Catch observation error (\code{logSDcatch}).
#' @param Fbarage Integer vector of length 2. Age range (inclusive) for reporting
#'   mean F (e.g., \code{c(2, maxage - 3)}).
#' @param nsurvey Integer. Number of surveys.
#' @param q Numeric length \code{nsurvey}. Survey catchability scalars.
#' @param surveyStart Numeric length \code{nsurvey}. Start of survey in seasonal
#'   time within a year (0–1).
#' @param surveyEnd Numeric length \code{nsurvey}. End of survey in seasonal
#'   time within a year (0–1).
#' @param surveySD Numeric length \code{nsurvey}. Observation sd on survey index
#'   (log scale).
#' @param SDR Numeric. Recruitment deviation sd (used with \code{recruitment.type}).
#' @param rhoR Numeric in [-1, 1]. AR(1) coefficient for recruitment deviations.
#'   If 0, recruitment deviations are i.i.d. normal.
#' @param recruitment.type Character. Deviation process type for recruitment:
#'   \code{"AR"} (default) or \code{"random"}.
#' @param F0 Numeric or length-\code{nyear} vector. Baseline fishing mortality:
#'   a scalar (expanded) for \code{"constant"} or the starting level for \code{"AR"}.
#' @param fishing.type Character. Fishing-mortality dynamics: \code{"constant"} or
#'   \code{"AR"}. When \code{"AR"}, an AR(1) process with sd \code{SDF} and
#'   autocorrelation \code{rho_F} is applied on the log scale.
#' @param SDF Numeric. Process sd for the fishing-mortality AR(1) on the log scale.
#' @param rho_F Numeric in [-1, 1]. AR(1) coefficient for fishing mortality.
#' @param SDM Numeric. Process sd for natural mortality dynamics
#'   (\code{mortality != "constant"}).
#' @param rho Numeric in [-1, 1]. AR(1) coefficient for natural mortality when
#'   \code{mortality == "AR"}.
#' @param propM Numeric scalar or array \code{[age \u00D7 year \u00D7 season]}.
#'   Fractional timing of natural mortality within the year/season (passed through).
#' @param propF Numeric scalar or array \code{[age \u00D7 year \u00D7 season]}.
#'   Fractional timing of fishing mortality within the year/season (passed through).
#' @param M_limit Numeric or \code{NA}. Upper (or lower for decreasing case) bound
#'   for time-varying M. Defaults to \code{2 * M} if \code{NA}.
#' @param mortality Character. Natural-mortality dynamics: \code{"constant"},
#'   \code{"rw"} (log-random walk with cap \code{M_limit}), \code{"AR"} (log-AR(1)),
#'   \code{"increase"} (linear +\code{SDM} per year up to \code{M_limit}), or
#'   \code{"decrease"} (linear -\code{SDM} per year down to \code{M_limit}, floored at 0).
#' @param recruitment Character. Stock–recruit relationship: \code{"BH_steep"}
#'   (Beverton–Holt with steepness), \code{"BH_env"} (BH with environmental effect
#'   \code{beta_env} and \code{env}), or other options the simulator accepts;
#'   value is passed through to the OM list.
#' @param env \code{NULL}, vector, or matrix. Environmental covariate(s) used when
#'   \code{recruitment == "BH_env"}. Will be coerced to a matrix if provided.
#' @param beta_env Numeric or \code{NULL}. Coefficient(s) for the environmental
#'   effect in \code{BH_env}. Required when \code{recruitment == "BH_env"}.
#' @param Fminage Integer. Minimum age fished (selectivity forced to 0 below this).
#' @param Fmaxage Integer. Maximum age fished (selectivity forced to 0 above this).
#' @param Qminage Integer vector length \code{nsurvey}. Minimum age selected by each survey.
#' @param Qmaxage Integer vector length \code{nsurvey}. Maximum age selected by each survey.
#' @param seed Integer or \code{NA}. RNG seed. If \code{NA}, a random seed is drawn.
#'
#' @details
#' \strong{Growth & weights:}
#' If \code{Linf} is numeric, von Bertalanffy lengths are generated over
#' \code{season_age = seq(0, maxage + 1/nseason, by = 1/nseason)} and converted
#' to weights using \code{w = a * L^b_size}. Otherwise, provide \code{weca}/\code{west}.
#'
#' \strong{Maturity & selectivity:}
#' Maturity and fishing selectivity follow logistic curves with slopes
#' \code{gamma}, \code{gamma_sel} and inflection ages \code{tau}, \code{tau_sel}.
#' Recruits (age 0) are forced to be immature.
#'
#' \strong{Fishing mortality (AR option):}
#' When \code{fishing.type = "AR"}, an AR(1) process is applied on the log scale:
#' \deqn{\omega_t = \rho_F \omega_{t-1} + \sqrt{1-\rho_F^2}\,\varepsilon_t,\quad
#' \varepsilon_t \sim N(0, SDF^2),\quad F_t = F_1 \exp(\omega_t - 0.5\,SDF^2).}
#'
#' \strong{Natural mortality:}
#' See \code{mortality} for options. For \code{"rw"} and \code{"AR"}, dynamics are
#' on the log scale; \code{"increase"}/\code{"decrease"} are linear on the natural
#' scale and truncated by \code{M_limit} and zero as appropriate.
#'
#' \strong{Recruitment deviations:}
#' If \code{recruitment.type = "AR"}, log-deviations follow AR(1) with sd \code{SDR}
#' and coefficient \code{rhoR}; if \code{rhoR = 0}, they become i.i.d. normal.
#'
#' \strong{Initial state:}
#' If \code{Fpast > 0}, an equilibrium-like \code{Ninit} is computed given
#' \code{M0[1]}, \code{Fpast}, and \code{sel}; otherwise \code{Ninit = NULL}.
#'
#' @return A named list of operating-model inputs. Important elements include:
#' \itemize{
#'   \item \code{west}, \code{weca}: weights for SSB and catch
#'   \item \code{mat}, \code{Fsel}, \code{Q}: maturity, fishing selectivity, survey selectivity
#'   \item \code{M}: time series of natural mortality
#'   \item \code{F0}: matrix of fishing mortality
#'   \item \code{years}, \code{age}, \code{nage}, \code{nyear}, \code{nseason}, \code{tEnd}
#'   \item Observation/variance parameters: \code{logSDcatch}, \code{logSDR}, \code{surveySD}
#'   \item Process controls and metadata: \code{rho}, \code{rhoR}, \code{SDF}, \code{SDM},
#'         \code{mortality}, \code{recruitment}, \code{h}, \code{env}, \code{beta_env}
#' }
#'
#' @seealso Functions in \pkg{smsR} that simulate or fit the OM using this list.
#'
#' @examples
#' # Minimal annual model with constant F and constant M
#' om <- sim_OM_parameters(
#'   nseason = 1, nyear = 30, Linf = 30, K = 0.6,
#'   M = 0.2, F0 = 0.25, fishing.type = "constant",
#'   mortality = "constant", R0 = 1e5, h = 0.8, seed = 42
#' )
#'
#' # Seasonal model with AR(1) F and AR(1) M
#' om_seas <- sim_OM_parameters(
#'   nseason = 4, nyear = 40, Linf = 35, K = 0.8,
#'   fishing.type = "AR", F0 = 0.2, SDF = 0.1, rho_F = 0.7,
#'   mortality = "AR", M = 0.25, SDM = 0.15, rho = 0.6,
#'   recruitment.type = "AR", SDR = 0.6, rhoR = 0.5,
#'   nsurvey = 2, q = c(0.6, 0.8), seed = 123
#' )
#'
#' @export
sim_OM_parameters <- function(nseason = 1,
                              nspace = 1,
                              nyear = 50,
                              F0 = 0.2,
                              Linf = 30,
                              weca = NULL,
                              west = NULL,
                              K = 1,
                              Fpast = 0,
                              t0 = -0.1,
                              gamma = 1.1,
                              tau = 5,
                              gamma_sel = 1.1,
                              tau_sel = 1,
                              maxage = 10,
                              a = 0.01,
                              b_size = 3,
                              M = 0.2,
                              b = 1,
                              R0 = 1e4,
                              theta = 2,
                              h = 0.8,
                              SDcatch = 0.01,
                              Fbarage = c(2, maxage-3),
                              nsurvey = 1,
                              q = rep(1, nsurvey),
                              surveySeason = rep(1, nsurvey),
                              surveyStart = rep(0, nsurvey),
                              surveyEnd = rep(0.99, nsurvey),
                              surveySD = rep(0.3, nsurvey),
                              SDR = .6,
                              rhoR = 0.7,
                              recruitment.type = 'AR',
                              rseason = 1,
                              fishing.type = 'AR',
                              SDF = 0.05,
                              rho_F = 0.7,
                              SDM = 0.1,
                              rho = 0.8,
                              propM = 0,
                              propF = 0,
                              M_limit = NA,
                              mortality = 'constant',
                              recruitment = 'BH_steep',
                              env = NULL,
                              beta_env = NULL,
                              Fminage = 0,
                              Fmaxage = maxage,
                              Qminage = rep(0, nsurvey),
                              Qmaxage = rep(maxage, nsurvey),
                              seed = NA
                              ){

  if(is.na(seed)){
    set.seed(round(runif(1, min = 1, max = 1e6)))
  }else{
    set.seed(seed)
  }

  if(is.na(M_limit)){
    M_limit <- 2*M
  }

  years <- 1:nyear
  nyear <- length(years)
  tEnd <- length(years)*nseason
  age <- 0:maxage

  ## Age stuff
  nage <- length(age)
  msel <- rep(1,nage)

  # Run growth

  if(is.numeric(Linf)){

  season_age <- seq(0, maxage+1/nseason, by = 1/nseason)

  Laa <- Linf*(1-exp(-K*(season_age-t0)))

  mat <- 1/(1+exp(-gamma*(season_age-tau)))
  mat[1] <- 0 # Don't let the recruits contribute to SSB in this model


  wage_ssb <- wage_catch <- wage_survey <- wage_mid <- array(NA, dim = c(nyear, nage, nseason))

  # Assign to matrices
  if(nseason == 1){

  wage_ssb <- t(replicate(nyear,a*Laa^b_size ))
  wage_catch <- t(replicate(nyear,a*Laa^b_size ))
  wage_survey <- t(replicate(nyear,a*Laa^b_size ))
  wage_mid <- t(replicate(nyear,a*Laa^b_size ))
  }else{
  sizes <- matrix(a*Laa^b_size, nrow = nage, ncol = nseason, byrow = TRUE)
  # Weight at age
  wage_tmp  <- replicate(nyear+1, sizes, simplify = 'array')
  wage_ssb <- wage_catch <- wage_survey <- wage_mid <- aperm(wage_tmp, c(1,3,2))

  # Maturity at age
  mat <- aperm(replicate(nyear + 1,matrix(mat, nrow = nage, ncol = nseason, byrow = TRUE),
                   simplify = 'array'), c(1,3,2))

    }
  }else{
  wage_ssb <- weca
  wage_catch <- west
  }


  sel <- 1/(1+exp(-gamma_sel*(age-tau_sel)))
  sel[age < Fminage] <- 0
  sel[age > Fmaxage] <- 0

  # Survey selectvitiy

  Q <- matrix(0, nage, nsurvey)

  if(length(q) < nsurvey){
    q <- rep(q, nsurvey)
  }

  for(i in 1:nsurvey){
  Q[,i] <- sel * q[i]
  Q[age < Qminage[i],i] <- 0
  Q[age > Qmaxage[i],i] <- 0
  }


  if(fishing.type == 'constant'){
    if(length(F0) == 1){
      F0 <- matrix(F0, nyear, nseason)
    }

  }


  # Do some magic for the seasonal AR models


  if(fishing.type == 'AR'){


    Fin <- rep(0, nyear*nseason)
    Fin[1] <- F0
    omega <- rep(NA, nyear*nseason)
    omega[1] <- 0

    for(i in 2:(nyear*nseason)){

      omega[i] <- rho_F*omega[i-1]+sqrt(1-rho_F^2)*rnorm(1,0,sd = sqrt(SDF))
      Fin[i]<- Fin[1]*exp(omega[i]-0.5*(SDF^2))


    }
    F0 <- matrix(Fin, nyear, nseason)
  }



  if(length(F0) != nyear*nseason){
    stop('wrong number of fishing mortality years')
  }


  if(fishing.type == 'est'){


  }

 #
  if(is.null(env) == FALSE){
    if(is.null(dim(env))){
      env <- as.matrix(env)
    }
  }

  # Option for natural mortality
  if(mortality == 'constant' & length(M) == 1){
    M0 <-    rep(M, nyear)
  }

  if(mortality == 'rw'){ # Natural mortality random  walk (this needs a bias adjuster)
    M0 <- rep(0, nyear)
    M0[1] <- M

    for(i in 2:nyear){
    M0[i] <- M0[i-1]*exp(rnorm(n = 1, mean = 0, sd = SDM))

    if(M0[i] > M_limit){ # Specify a max limit for natural mortliaty
      M0[i] <- M_limit
    }

    }


  }

  if(mortality == 'AR'){ # Natural mortality autocorrelation

    M0 <- rep(0, nyear)
    M0[1] <- M
    omega <- rep(NA, nyear)
    omega[1] <- 0

    for(i in 2:nyear){

        omega[i] <- rho*omega[i-1]+sqrt(1-rho^2)*rnorm(1,0,sd = sqrt(SDM))
        M0[i]<-M0[1]*exp(omega[i]-0.5*(SDM^2))


    }

  }




  if(mortality == 'increase'){
    M0 <- rep(0, nyear)
    M0[1] <- M

    for(i in 2:nyear){
      M0[i] <- M0[i-1]+SDM

      if(M0[i] > M_limit){ # Specify a max limit for natural mortliaty
        M0[i] <- M_limit
      }

    }

  }


  if(mortality == 'decrease'){

    M0 <- rep(0, nyear)
    M0[1] <- M

    for(i in 2:nyear){
      M0[i] <- M0[i-1]-SDM

      if(M0[i] < M_limit){ # Specify a max limit for natural mortliaty
        M0[i] <- M_limit
      }
      M0[M0 < 0] <- 0

    }
  }


  if(length(M0) != nyear){
    stop('wrong number of natural mortality years')
  }



  # Insert Rdev function here
  Rdev <- NA

  if(rhoR == 0){
    recruitment.type <- 'random'
  }


  if(recruitment.type  == 'random'){
    Rdev <- rnorm(nyear, mean = 0, sd = SDR)
  }

  if(recruitment.type == 'AR'){

    Rdev <- rep(0, nyear)
    Rdev[1] <- 0#rnorm(1, mean = 0, sd = SDR)

    omegaR <- rep(NA, nyear)
    omegaR[1] <- 0

    for(i in 2:nyear){

      omegaR[i] <- rhoR*omegaR[i-1]+sqrt(1-rhoR^2)*rnorm(1,mean = 0,sd = SDR)
      #Rdev[i]<- omegaR[i]-0.5*(SDR^2)

    #print('test')
    }
    Rdev <- omegaR


  }



  if(length(b) == 1){
  bout <- rep(b, nyear)
  }



  # if(recruitment == 'BH'){
  #   R0 <- alpha/beta
  # }
  #
  # if(recruitment == 'Ricker'){
  #   R0 <- 1/beta
  # }
#
#   Custom initial distribution
  if(Fpast > 0){

  Ninit <- matrix(0, nage)

  Ninit[1:(nage - 1)] <- R0 * exp(-(age[1:(nage - 1)] * M0[1]+sel[1:(nage - 1)] * age[1:(nage - 1)] * Fpast )) # Only valid for equal M2 per age
  Ninit[nage] <- R0 * exp(-(M0[1] * age[nage]+sel[nage] * Fpast * age[nage])) / (1 - exp(-(M0[1]+Fpast*sel[nage])))
  }else{
    Ninit = NULL
  }



  if(nseason == 1){

    Fnseason <- matrix(rep(1, nspace))

  }else{

    if(exists('Fnseason') !=0)
    warning('Fisheries per season not provided. Assuming equal')
    Fnseason <- matrix(rep(1, nspace, nseason))/nseason

  }


  if(exists('omega') == FALSE){
    omega = NA
  }


  if(length(propM) == 1){
    propM = array(propM, dim = c(nage, nyear, nseason))
  }
  if(length(propF) == 1){
    propF = array(propF, dim = c(nage, nyear, nseason))
  }

  if(recruitment == 'BH_steep' | recruitment == 'BH_env'){
    h <- h
  }else{
    h <- NULL
  }

  if(recruitment == 'BH_env'){
    if(is.null(beta_env)) stop('Provide environmental coefficient beta_env')
  }



  df <-list(      #### Parameters #####
                  west = wage_ssb,
                  weca = wage_catch,
                  Msel = msel,
                  mat = mat,
                  nage = nage,
                  age = age,
                  tau = tau,
                  nseason = nseason,
                  nyear = nyear,
                  tEnd = tEnd,
                  logq = log(q),
                  surveySeason = surveySeason,
                  surveyStart = surveyStart,
                  surveyEnd = surveyEnd,
                  M = M0,
                  K = K,
                  Ninit = Ninit,
                  Linf = Linf,
                  nsurvey = nsurvey,
                  Q = Q,# Frequency of survey years (e.g., 2 is every second year)
                  # variance parameters
                  logSDcatch = log(SDcatch),
                  logSDR = log(SDR), # Fixed in stock assessment ,
                  surveySD = surveySD,
                  years = years,
                  b = bout,
                  h = h,
                  nspace = nspace,
                  Fmodel = 'sim', # Fix later maybe
                  move = FALSE,
                  Fbarage = Fbarage,
                  propM = propM,
                  propF = propF,
                  Fsel = sel,
                  Fnseason = Fnseason,
                  F0 = F0,
                  rseason = rseason,
                  theta = theta,
                  SDM = SDM,
                  rho = rhoR,
                  omega = omega,
                  M_limit = M_limit,
                  mortality = mortality,
                  recruitment = recruitment,
                  env = env,
                  beta_env = beta_env,
                  R0 = R0,
                  Fminage = Fminage,
                  Fmaxage = Fmaxage,
                  Qminage = Qminage,
                  Qmaxage = Qmaxage



                  # Parameters from the estimation model
  )


  return(df)

}
