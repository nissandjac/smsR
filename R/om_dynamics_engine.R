#' Shared population-dynamics engine for smsR operating models
#'
#' Internal helper that runs the season/year loop (recruitment, mortality,
#' movement, catch and survey calculation) shared by \code{\link{run.agebased.sms.op}}
#' (historical OM run) and \code{\link{forecast_op}} (OM forecast). Callers are
#' responsible for building the time-indexed input arrays and for supplying an
#' \code{Fseason_fun} closure that encodes how fishing mortality is derived in
#' their context (extracted from \code{df} for the historical run, user-supplied
#' for the forecast).
#'
#' @param df List of OM parameters (recruitment settings, proportions, survey
#'   config, etc. - the season loop reads these directly from \code{df}).
#' @param nyear Integer. Number of years to simulate.
#' @param year,year_1 Year labels (length \code{nyear} and \code{nyear + 1}) used as dimnames.
#' @param nage,age,maxage Age dimension, age labels, and maximum age.
#' @param nseason,nspace Number of seasons and spatial areas.
#' @param M,weca,west,mat,selin Time-indexed arrays, dim \code{c(nage, >=nyear, 1, nseason)}.
#' @param movemat Movement array, dim \code{c(nage, nyear, 1, nseason)}.
#' @param Ninit Numbers-at-age used to initialize year 1.
#' @param R0 Unfished recruitment, used by the stock-recruit formulas.
#' @param SSB_0 Unfished SSB (NA if not used by the chosen recruitment model).
#' @param Fseason_fun function(yr, season, space, sel_season) returning the
#'   realized F-at-age vector for that season/space.
#' @param Rin_fun function(yr, space) returning recruitment when
#'   \code{df$recruitment == "estimated"} (unused for other recruitment models).
#'   The historical run indexes \code{df$Rin} directly by \code{yr}; a forecast
#'   cannot do that (there is no future \code{df$Rin}) and must supply its own
#'   closure, e.g. a geometric mean of recent historical recruitment.
#' @param env Optional environmental covariate matrix, indexed as \code{env[yr, ]}
#'   (already aligned to the engine's internal \code{yr} index by the caller).
#' @param stochastic Logical. If \code{FALSE}, recruitment deviations are drawn as zero.
#'
#' @return A list of the raw (noise-free) simulation outputs shared by both callers.
#' @keywords internal
simulate_om_dynamics <- function(df,
                                  nyear,
                                  year,
                                  year_1,
                                  nage,
                                  age,
                                  maxage,
                                  nseason,
                                  nspace,
                                  M,
                                  weca,
                                  west,
                                  mat,
                                  selin,
                                  movemat,
                                  Ninit,
                                  R0,
                                  SSB_0,
                                  Fseason_fun,
                                  Rin_fun = NULL,
                                  env = NULL,
                                  stochastic = TRUE) {

  SDR <- exp(df$logSDR)

  draw_err <- function(sd) {
    if (stochastic) rnorm(1, mean = 0, sd = sd) else 0
  }

  SSB <- matrix(NA, nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  SSB.all <- array(NA, dim = c(nyear, nspace, nseason), dimnames = list(year = year, space = 1:nspace, season = 1:nseason))
  SSB.weight <- matrix(NA, nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  Catch <- matrix(NA, nyear, dimnames = list(year = year))
  Catch.age <- matrix(NA, nage, nyear, dimnames = list(age = age, year = year))
  CatchN <- matrix(NA, nyear, dimnames = list(year = year))
  CatchN.age <- matrix(NA, nage, nyear, dimnames = list(age = age, year = year))

  R.save <- matrix(NA, nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  R.err.save <- matrix(NA, nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  Fseason.save <- array(NA, dim = c(nage, nyear, nspace, nseason), dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))

  N.save.age <- array(0, dim = c(nage, nyear + 1, nspace, nseason), dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  N.save.age.mid <- array(NA, dim = c(nage, nyear + 1, nspace, nseason), dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  V.save <- array(NA, dim = c(nyear, nspace, nseason), dimnames = list(year = year, space = 1:nspace, season = 1:nseason))
  E.save <- array(NA, dim = c(nyear, nspace, nseason), dimnames = list(year = year, space = 1:nspace, season = 1:nseason))

  Catch.save.age <- array(0, dim = c(nage, nyear, nspace, nseason), dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))
  CatchN.save.age <- array(0, dim = c(nage, nyear, nspace, nseason), dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))

  survey <- array(NA, dim = c(nage, nyear, df$nsurvey))
  survey.true <- array(NA, dim = c(nage, nyear, nspace, df$nsurvey), dimnames = list(age = 0:maxage, year = year, space = 1:nspace, survey = 1:df$nsurvey))

  age_comps_catch_space <- array(NA, dim = c(maxage, nyear, nspace), dimnames = list(age = 1:maxage, year = year, space = 1:nspace))
  age_comps_OM <- array(NA, dim = c(nage, nyear, nspace, nseason), dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))

  Z.save <- array(NA, dim = c(nage, nyear, nspace, nseason), dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))

  for (space in 1:nspace) {
    N.save.age[, 1, space, 1] <- Ninit / nspace
    N.save.age.mid[, 1, space, 1] <- N.save.age[, 1, space, 1] * exp(-0.5 * (M[, 1, space, 1] / nseason))
  }

  for (yr in 1:nyear) {

    w_catch <- weca[, yr, , , drop = FALSE]
    w_catch[is.na(w_catch)] <- 0

    w_ssb <- west[, yr, , , drop = FALSE]
    w_ssb[is.na(w_ssb)] <- 0

    sel <- selin[, yr, , , drop = FALSE]
    sel[is.na(sel)] <- 0

    Myear <- M[, yr, , , drop = FALSE]
    mat.year <- mat[, yr, , , drop = FALSE]

    for (season in 1:nseason) {
      for (space in 1:nspace) {

        Fseason <- Fseason_fun(yr, season, space, sel[, 1, space, season])

        if (season == 1) {
          SSB.weight[yr, space] <- sum(N.save.age[, yr, space, 1] * w_ssb[, , space, season] *
            mat.year[, 1, space, season] * exp(-(Myear[, 1, space, season] * df$propM[, yr, season] + Fseason * df$propF[, yr, season])), na.rm = TRUE)
          SSB[yr, space] <- SSB.weight[yr, space]

          SSB.all[1, space, 1] <- sum(N.save.age[, 1, space, 1] * mat.year[, 1, space, season] *
            exp(-(Myear[, 1, space, season] * df$propM[, yr, season] + Fseason * df$propF[, yr, season])), na.rm = TRUE)
        }

        if (season == df$rseason) { # Recruitment season

          if (df$recruitment == "Ricker") {
            R <- df$alpha * SSB[yr, ] * exp(-betaSR * SSB[yr, ])
          }

          if (df$recruitment == "hockey") {
            R <- log(df$R0 / df$betaSR) + log(SSB[yr, space])

            if (SSB[yr, space] >= df$betaSR) {
              R <- log(df$R0)
            }

            R <- exp(R)

            Ry <- draw_err(SDR)
            R.err <- exp(-0.5 * df$b[yr] * SDR^2 + Ry)

            N.save.age[1, yr, space, season] <- R * R.err
            R.save[yr, space] <- R * R.err
            R.err.save[yr, space] <- R.err
          }

          if (df$recruitment == "Rmean") {

            R <- df$R0

            Ry <- draw_err(SDR)
            R.err <- exp(-0.5 * df$b[yr] * SDR^2 + Ry)

            N.save.age[1, yr, space, season] <- R * R.err
            R.save[yr, space] <- R * R.err
          }

          if (df$recruitment == "estimated") {
            R <- Rin_fun(yr, space)

            N.save.age[1, yr, space, season] <- R
            R.save[yr, space] <- R
          }

          if (df$recruitment == "BH_steep") {
            err <- draw_err(SDR)

            R <- (4 * df$h * R0 * SSB[yr, space] / (SSB_0 * (1 - df$h) + SSB[yr, space] * (5 * df$h - 1))) *
              exp(-0.5 * SDR^2 + err) # Should probably add recruitment bias correction

            N.save.age[1, yr, space, season] <- R
            R.save[yr, space] <- R
            R.err.save[yr, space] <- err
          }

          if (df$recruitment == "BH_env") {
            err <- draw_err(SDR)

            env_tot <- sum(df$beta_env * env[yr, ])

            R <- (4 * df$h * R0 * SSB[yr, space] / (SSB_0 * (1 - df$h) + SSB[yr, space] * (5 * df$h - 1))) *
              exp(-0.5 * SDR^2 + err + env_tot)

            N.save.age[1, yr, space, season] <- R
            R.save[yr, space] <- R
            R.err.save[yr, space] <- err
          }
        }

        Mseason <- Myear[, 1, space, season]

        Z <- Mseason + Fseason

        if (Z[1] == 0) {
          Z[1] <- Z[2] # To avoid NaN calc
        }

        Z.save[, yr, space, season] <- Z
        Fseason.save[, yr, space, season] <- Fseason

        if (((space - 1) == 0)) {
          spaceidx <- 2
        }
        if (space == nspace) {
          spaceidx <- nspace - 1
        }
        if (space > 1 & space < nspace) {
          spaceidx <- c(space - 1, space + 1)
        }

        if (df$move == FALSE) {
          spaceidx <- 1
        }

        if (season < nseason) {

          for (k in 1:length(spaceidx)) {

            if (spaceidx[k] > 1 & (spaceidx[k] < nspace)) {
              space_multiplier <- 0.5
            } else {
              space_multiplier <- 1
            }
            Nin.tmp <- N.save.age[, yr, spaceidx[k], season] * exp(-Z) * (movemat[, yr, spaceidx[k], season]) * space_multiplier

            if (k == 1) {
              Nin <- Nin.tmp
            } else {
              Nin <- Nin + Nin.tmp
            }
          }

          N.save.age[, yr, space, season + 1] <- N.save.age[, yr, space, season] * exp(-Z) -
            N.save.age[, yr, space, season] * exp(-Z) * (movemat[, yr, space, season]) +
            Nin

          age_comps_OM[, yr, space, season] <- N.save.age[, yr, space, season] / sum(N.save.age[, yr, space, season])

          SSB.all[yr, space, season] <- sum(N.save.age[, yr, space, season] * mat.year[, 1, space, season], na.rm = T)
          V.save[yr, space, season] <- sum(N.save.age[, yr, space, season] * sel[, 1, space, season] * w_catch[, 1, space, season])

          if (max(Fseason) > 0) {
            Zcatch <- Z
            Zcatch[1] <- Zcatch[2]

            Catch.save.age[, yr, space, season] <- (Fseason / (Zcatch)) * (1 - exp(-(Zcatch))) * N.save.age[, yr, space, season] * w_catch[, , space, season]
            CatchN.save.age[, yr, space, season] <- (Fseason / (Zcatch)) * (1 - exp(-(Zcatch))) * N.save.age[, yr, space, season]

            E.save[yr, space, season] <- sum(Catch.save.age[, yr, space, season]) / V.save[yr, space, season]
          }
        } else {
          for (k in 1:length(spaceidx)) {
            if ((spaceidx[k]) > 1 & (spaceidx[k] < nspace)) {
              space_multiplier <- 0.5
            } else {
              space_multiplier <- 1
            }

            Nin.tmp <- N.save.age[1:(nage - 2), yr, spaceidx[k], season] * exp(-Z[1:(nage - 2)]) * (movemat[1:(nage - 2), yr, spaceidx[k], season]) * space_multiplier

            Nin.plus.tmp <- (N.save.age[nage - 1, yr, spaceidx[k], nseason] * exp(-Z[nage - 1]) +
              N.save.age[nage, yr, spaceidx[k], nseason] * exp(-Z[nage])) *
              (movemat[nage, yr, spaceidx[k], season]) * space_multiplier

            if (k == 1) {
              Nin <- Nin.tmp
              Nin.plus <- Nin.plus.tmp
            } else {
              Nin <- Nin + Nin.tmp
              Nin.plus <- Nin.plus.tmp + Nin.plus
            }
          }

          N.save.age[2:(nage - 1), yr + 1, space, 1] <- N.save.age[1:(nage - 2), yr, space, season] * exp(-Z[1:(nage - 2)]) -
            N.save.age[1:(nage - 2), yr, space, season] * exp(-Z[1:(nage - 2)]) * (movemat[1:(nage - 2), yr, space, season]) +
            Nin
          # Plus group
          Nsurvive.plus <- (N.save.age[nage - 1, yr, space, nseason] * exp(-Z[nage - 1]) +
            N.save.age[nage, yr, space, nseason] * exp(-Z[nage]))

          Nout.plus <- Nsurvive.plus * (movemat[nage, yr, space, season])

          N.save.age[nage, yr + 1, space, 1] <- Nsurvive.plus - Nout.plus + Nin.plus

          age_comps_OM[, yr, space, season] <- N.save.age[, yr, space, season] / sum(N.save.age[, yr, space, season])

          SSB.all[yr, space, season] <- sum(N.save.age[, yr, space, season] * mat.year[, 1, space, season], na.rm = T)
          V.save[yr, space, season] <- sum(N.save.age[, yr, space, season] * sel[, 1, space, season] * w_catch[, 1, space, season])

          if (max(Fseason) > 0) {
            Zcatch <- Z
            Zcatch[1] <- Zcatch[2]

            Catch.save.age[, yr, space, season] <- (Fseason / (Zcatch)) * (1 - exp(-(Zcatch))) * N.save.age[, yr, space, season] * w_catch[, 1, space, season]
            CatchN.save.age[, yr, space, season] <- (Fseason / (Zcatch)) * (1 - exp(-(Zcatch))) * N.save.age[, yr, space, season]

            E.save[yr, space, season] <- sum(Catch.save.age[, yr, space, season]) / V.save[yr, space, season]
          }
        }
      }
    } # End of season loop

    if (nseason > 1) {
      Catch.age[, yr] <- apply(Catch.save.age[, yr, , ], MARGIN = 1, FUN = sum)
      Catch[yr] <- sum(Catch.save.age[, yr, , ])

      CatchN.age[, yr] <- apply(CatchN.save.age[, yr, , ], MARGIN = 1, FUN = sum)
      CatchN[yr] <- sum(CatchN.save.age[, yr, , ])
    } else {
      if (nspace == 1) {
        Catch.age[, yr] <- Catch.save.age[, yr, , ]
        Catch[yr] <- sum(Catch.save.age[, yr, , ])

        CatchN.age[, yr] <- CatchN.save.age[, yr, , ]
        CatchN[yr] <- sum(CatchN.save.age[, yr, , ])
      } else {
        Catch.age[, yr] <- rowSums(Catch.save.age[, yr, , ])
        Catch[yr] <- sum(Catch.save.age[, yr, , ])

        CatchN.age[, yr] <- rowSums(CatchN.save.age[, yr, , ])
        CatchN[yr] <- sum(CatchN.save.age[, yr, , ])
      }
    }

    for (surv in 1:df$nsurvey) {
      if (df$move == FALSE) {
        if (df$surveyEnd[surv] == 0) {
          survey[, yr, surv] <- N.save.age[, yr, , df$surveySeason[surv]] *
            exp(-Z.save[, yr, space, df$surveySeason[surv]]) * df$Q[, surv]
        } else {
          Ntmp.s <- N.save.age[, yr, , df$surveySeason[surv]] * (exp(-Z.save[, yr, space, df$surveySeason[surv]] * df$surveyStart[surv]))
          survey[, yr, surv] <- Ntmp.s * (1 - exp(-Z.save[, yr, space, df$surveySeason[surv]] * (df$surveyEnd[surv] - df$surveyStart[surv]))) /
            (Z.save[, yr, space, df$surveySeason[surv]] * (df$surveyEnd[surv] - df$surveyStart[surv])) * df$Q[, surv]
        }
        # move == FALSE collapses to a single effective space; mirror the
        # aggregate survey into survey.true so it isn't left all-NA.
        survey.true[, yr, , surv] <- survey[, yr, surv]
      } else {
        for (space in 1:nspace) {
          if (df$surveyEnd[surv] == 0) {
            survey.true[, yr, space, surv] <- N.save.age[, yr, space, df$surveySeason[surv]] *
              exp(-Z.save[, yr, space, df$surveySeason[surv]]) * df$Q[, surv]
          } else {
            Ntmp.s <- N.save.age[, yr, space, df$surveySeason[surv]] * (exp(-Z.save[, yr, space, df$surveySeason[surv]] * df$surveyStart[surv]))
            survey.true[, yr, space, surv] <- Ntmp.s * (1 - exp(-Z.save[, yr, space, df$surveySeason[surv]] * (df$surveyEnd[surv] - df$surveyStart[surv]))) /
              (Z.save[, yr, space, df$surveySeason[surv]] * (df$surveyEnd[surv] - df$surveyStart[surv])) * df$Q[, surv]
          }
        }
      }

      survey[survey == 0] <- -1 # For TMB
    }

    for (space in 1:nspace) {
      if (nseason > 1) {
        Catch.tmp <- rowSums(CatchN.save.age[, yr, space, ])
      } else {
        Catch.tmp <- CatchN.save.age[, yr, space, ]
      }

      Catch.tot <- sum(CatchN.save.age[, yr, space, ])

      age_comps_catch_space[1:(maxage - 1), yr, space] <- Catch.tmp[2:(maxage)] / Catch.tot
      age_comps_catch_space[maxage, yr, space] <- sum(Catch.tmp[(maxage + 1):nage]) / Catch.tot
    }
  } # End of year loop

  if (df$move == FALSE) {
    Nsave <- N.save.age[, , , nspace]
  } else {
    Nsave <- apply(N.save.age[, , , 1], 2, rowSums)
  }

  if (df$move == TRUE) {
    for (i in 1:df$nsurvey) {
      survey[, , i] <- apply(survey.true[, , , i], MARGIN = c(1, 2), sum)
    }
  }

  Fbar <- rep(0, nyear)

  for (time in 1:nyear) {
    Fbar[time] <- sum(Fseason.save[(df$Fbarage[1] + 1):(df$Fbarage[2] + 1), time, , ]) / (df$Fbarage[2] - df$Fbarage[1] + 1)
  }

  list(
    N.save = Nsave,
    SSB = SSB,
    N.save.age = N.save.age,
    R.save = R.save,
    R.err.save = R.err.save,
    V.save = V.save,
    E.save = E.save,
    SSB.all = SSB.all,
    Catch.save.age = Catch.save.age,
    CatchN.save.age = CatchN.save.age,
    Catch = Catch,
    Catch.age = Catch.age,
    CatchN = CatchN,
    CatchN.age = CatchN.age,
    survey = survey,
    survey.true = survey.true,
    Fbar = Fbar,
    M = M,
    age_comps_OM = age_comps_OM,
    age_catch = age_comps_catch_space,
    Z = Z.save,
    Fseason = Fseason.save
  )
}
