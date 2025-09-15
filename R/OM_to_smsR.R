#' Turn an operating model into smsR-ready inputs
#'
#' @description
#' Helper that reshapes/sanitizes outputs from an operating model (**OM**) into
#' a compact list of arrays and options consumed by **smsR**. It builds quarter-
#' level catch and survey inputs and (optionally) perturbs catch with lognormal
#' noise using `df$logSDcatch`.
#'
#' @param OM list
#'   Operating model output. Must include at least:
#'   - `survey` : array or list of survey observations (dimensions must align with
#'     `nage × nyear × nseason`, or be coercible).
#'   - `CatchN.save.age` : 4D array of catch-at-age numbers with dimensions
#'     `[age, year, season, space]`; smsR uses the first spatial stratum
#'     (i.e. `[..., , , 1]`).
#'
#' @param df list
#'   Operating-model/meta parameters required for shaping inputs; must include:
#'   - `nyear` (integer): number of years.
#'   - `nage` (integer): number of ages.
#'   - `nseason` (integer): number of seasons (quarters).
#'   - `logSDcatch` (numeric, log-scale): catch observation SD used when adding
#'     lognormal noise (see Details). May be length-1 or length-`nage`.
#'
#' @param recmodel integer(1)
#'   Recruitment model code for **estimation** in smsR. Typical choices:
#'   - `0` = none/fixed,
#'   - `1` = Beverton–Holt,
#'   - `2` = Ricker.
#'   (Adjust if your package uses different codes.)
#'
#' @param MCV integer or NULL
#'   Optional grouping index for natural mortality (**M**) by age (e.g.,
#'   a vector of length `nage` mapping ages to M groups). Use `NULL` to
#'   treat ages independently.
#'
#' @param M_max numeric(1)
#'   Upper bound used when defining priors/penalties on natural mortality.
#'   Defaults to `max(df$age)` if available (adapt to your convention).
#'
#' @return list
#'   A named list ready for smsR, typically including:
#'   - `Catchobs` : `array[nage, nyear, nseason]` catch-at-age numbers.
#'   - `Catchobs_err` : same shape as `Catchobs`, with optional lognormal noise.
#'   - `Surveyobs` : survey observations aligned to `nage × nyear × nseason`.
#'   - `recmodel`, `MCV`, `M_max` : passed-through control parameters.
#'   Additional elements may be included as needed by smsR.
#'
#' @details
#' **Catch noise:** If `df$logSDcatch` is provided, `Catchobs_err` is formed as
#' `Catchobs * exp(ε)`, where `ε ~ N(0, sd = exp(df$logSDcatch))`. If
#' `logSDcatch` is length-`nage`, age-specific noise is applied per year/season.
#'
#' **Dimensionality:** `CatchN.save.age` is expected as
#' `[age = df$nage, year = df$nyear, season = df$nseason, space >= 1]`. The
#' function uses the first spatial dimension. `Surveyobs` must be conformable to
#' `nage × nyear × nseason`; otherwise, provide a coercible structure.
#'
#' @seealso
#' \code{\link{smsR}} model estimation functions; recruitment helpers for
#' Beverton–Holt/Ricker; data constructors for surveys and catches.
#'
#' @examples
#' # Minimal mock example (dimensions only)
#' set.seed(1)
#' df <- list(nyear = 3, nage = 5, nseason = 4, logSDcatch = log(0.2))
#' OM <- list(
#'   survey = array(runif(df$nage * df$nyear * df$nseason),
#'                  dim = c(df$nage, df$nyear, df$nseason)),
#'   CatchN.save.age = array(rpois(df$nage * df$nyear * df$nseason * 1, 100),
#'                           dim = c(df$nage, df$nyear, df$nseason, 1))
#' )
#' out <- OM_to_smsR(OM, df, recmodel = 1, MCV = NULL, M_max = df$nage)
#' str(out)
#'
#' @export
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
