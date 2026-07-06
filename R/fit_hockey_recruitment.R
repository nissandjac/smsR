#' Fit a hockey-stick stock-recruit relationship to an assessment's history
#'
#' Internal helper shared by \code{\link{getFmsy}} and \code{\link{get_OM_parameters}}
#' so both use the exact same fitted breakpoint/alpha when the underlying
#' assessment (\code{recmodel == 2}) has no analytic hockey-stick SR relationship
#' to draw on, and one has to be fit post-hoc to the assessment's own (SSB, R)
#' history instead.
#'
#' @param df.tmb list of smsR input data
#' @param sas fitted smsR model
#' @param method \code{"nls"} (default) or \code{"segmented"} - see \code{\link{getFmsy}}.
#'
#' @return list with \code{alphaSR}, \code{betaSR}, and the \code{(SSB, R)} pairs used to fit them.
#' @keywords internal
.fit_hockey_recruitment <- function(df.tmb, sas, method = "nls") {
  r_obs   <- getR(df.tmb, sas)$R
  ssb_obs <- getSSB(df.tmb, sas)$SSB[seq_along(r_obs)]

  if (method == "nls") {
    hs_nll <- function(pars) {
      alpha <- exp(pars[1])
      Bknee <- exp(pars[2])
      r_pred <- alpha * pmin(ssb_obs, Bknee)
      sum((log(r_obs) - log(r_pred))^2)
    }

    init_pars <- c(log(mean(r_obs) / mean(ssb_obs)), log(median(ssb_obs)))
    hs_opt <- optim(init_pars, hs_nll, method = "Nelder-Mead")

    alphaSR <- exp(hs_opt$par[1])
    betaSR  <- exp(hs_opt$par[2])
  } else if (method == "segmented") {
    if (!requireNamespace("segmented", quietly = TRUE)) {
      stop("Package 'segmented' is required for method = 'segmented'. Install it with install.packages('segmented').")
    }

    lin_fit <- lm(r_obs ~ ssb_obs)
    seg_fit <- segmented::segmented(lin_fit, seg.Z = ~ssb_obs, psi = median(ssb_obs))

    betaSR  <- as.numeric(seg_fit$psi[, "Est."])
    alphaSR <- as.numeric(segmented::slope(seg_fit)$ssb_obs["slope1", "Est."])
  } else {
    stop("method must be one of: nls, segmented")
  }

  list(alphaSR = alphaSR, betaSR = betaSR, r_obs = r_obs, ssb_obs = ssb_obs)
}
