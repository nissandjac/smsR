#' Calculate Fmsy from a sms model
#'
#' @param sas Fitted sms model
#' @param nyears number of years to simulate
#' @param recruitment recruitment type to calculate Fmsy. Currently only hockey works.
#' @param plotMSY TRUE OR FALSE to plot the result. When \code{recruitment == "hockey"} and
#'   \code{sas$dat$recmodel == 2}, this also plots the observed (SSB, R) pairs with the fitted
#'   hockey-stick curve overlaid, so the fit can be checked visually.
#' @param stochastic stochastic fmsy?
#' @param nruns number of replicates in stochastic fmsy
#' @param method Only used when \code{recruitment == "hockey"} and \code{sas$dat$recmodel == 2}
#'   (i.e. recruitment was estimated freely per year in the assessment, so the hockey-stick
#'   has to be fit post-hoc to the assessment's own (SSB, R) history). \code{"nls"} (default)
#'   fits alpha/beta by minimizing squared log-residuals with \code{optim()}. \code{"segmented"}
#'   instead uses the \pkg{segmented} package: an OLS line is fit and then refined into a
#'   breakpoint model via \code{segmented::segmented()}, using the fitted breakpoint as beta
#'   and the pre-breakpoint segment's slope as alpha. Requires the \pkg{segmented} package.
#'
#' @return returns the fishing mortality that leads to Fmsy
#' @export
#'
getFmsy <- function(sas,
                    nyears = 50,
                    recruitment = "hockey",
                    plotMSY = FALSE,
                    stochastic = FALSE,
                    nruns = 1000,
                    method = "nls") {
  #  Number
  df.tmb <- sas$dat
  years <- 1:nyears
  nyears <- length(years)

  # Get the index for forecast
  if(dim(df.tmb$weca)[2] == df.tmb$nyears){
    fc_idx <- df.tmb$nyears
  }else{
    fc_idx <- df.tmb$nyears+1
  }


  if(sas$dat$randomM == 1){

  M_estimated <- getM_var(sas$dat, sas) %>% filter(years %in% max(sas$dat$years))


  M <-   array(M_estimated$M_var, dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  }else{
  M <- array(df.tmb$M[, fc_idx, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  }

  mat_last <- df.tmb$Mat[, fc_idx, , drop = FALSE]   # nage × 1 × nseason
  mat <- mat_last[, rep(1, nyears), , drop = FALSE]


  weca_last <- df.tmb$weca[, fc_idx, , drop = FALSE]   # nage × 1 × nseason
  weca <- weca_last[, rep(1, nyears), , drop = FALSE]

  west_last <- df.tmb$west[, fc_idx, , drop = FALSE]   # nage × 1 × nseason
  west <- west_last[, rep(1, nyears), , drop = FALSE]

  propM_last <- df.tmb$propM[, fc_idx, , drop = FALSE]   # nage × 1 × nseason
  propM <- propM_last[, rep(1, nyears), , drop = FALSE]

  propF_last <- df.tmb$propF[, fc_idx, , drop = FALSE]   # nage × 1 × nseason
  propF <- propF_last[, rep(1, nyears), , drop = FALSE]

  Fselin <- getSelectivity(df.tmb, sas)
  #Fselin <- Fselin[Fselin$block == max(Fselin$block),]
  # I like this scaled to 1 for consistency

  # blocks <- unique(Fselin$blocks)


  # for(i in 1:length(blocks)){
  #
  #   Fsel$Fsel[Fsel$blocks == blocks[i]] <- Fsel$Fsel[Fsel$blocks == blocks[i]]/max(Fsel$Fsel[Fsel$blocks == blocks[i]])
  #
  # }

  for (i in 1:df.tmb$nseason) {
    if (i == 1) {
      Fsel <- array(Fselin$selec[Fselin$season == i], dim = c(df.tmb$nage, nyears, 1))
    } else {
      tmp <- array(Fselin$selec[Fselin$season == i], dim = c(df.tmb$nage, nyears, 1))
      Fsel <- abind::abind(Fsel, tmp, along = 3)
    }
  }

  parms.est <- getEstimatedParms(sas)

  Ninit <- exp(parms.est$value[parms.est$parameter == "logNinit"])
  # SDR <- exp(parms.est$value[parms.est$parameter == 'logSDrec'])

  if (recruitment == "hockey") {


    if(sas$dat$recmodel == 1){
    alphaSR <- exp(parms.est$value[parms.est$parameter == "logalpha"])
    betaSR <- df.tmb$betaSR
    R0 <- alphaSR * df.tmb$betaSR
    h_val <- NA
    Ninit <- c(R0,exp(parms.est$value[parms.est$parameter == 'logNinit']))
    }

    if(sas$dat$recmodel == 2){
      hs_fit  <- .fit_hockey_recruitment(sas$dat, sas, method = method)
      alphaSR <- hs_fit$alphaSR
      betaSR  <- hs_fit$betaSR

      if (plotMSY == TRUE) {
        ssb_seq <- seq(0, max(hs_fit$ssb_obs) * 1.05, length.out = 200)
        hs_fit_df <- data.frame(SSB = ssb_seq, R = alphaSR * pmin(ssb_seq, betaSR))
        hs_obs_df <- data.frame(SSB = hs_fit$ssb_obs, R = hs_fit$r_obs)

        print(
          ggplot() +
            geom_point(data = hs_obs_df, aes(x = SSB, y = R), alpha = 0.6) +
            geom_line(data = hs_fit_df, aes(x = SSB, y = R), color = "steelblue", linewidth = 1.1) +
            geom_vline(xintercept = betaSR, linetype = "dashed", color = "gray40") +
            theme_classic() +
            scale_x_continuous("SSB") +
            scale_y_continuous("Recruitment") +
            ggtitle(paste0("Hockey-stick fit (method = '", method, "')"))
        )
      }

      R0      <- alphaSR * betaSR
      h_val   <- NA
      Ninit   <- c(R0, exp(parms.est$value[parms.est$parameter == 'logNinit']))
    }

  } else {
    alphaSR <- NA
    betaSR <- NA
    R0 <- exp(parms.est$value[parms.est$parameter == "logR0"])
    h_val <- attr(sas$obj$env$parameters$logh, "shape")
    Ninit <- NULL
  }

  if(sas$dat$nenv >0){

    beta_env <- sas$reps$par.fixed[names(sas$reps$par.fixed) == 'env']

    # For environmental input, repeat the last year
    env <- matrix(df.tmb$env_matrix[df.tmb$nyears], ncol = df.tmb$nenv, nrow = nyears)

  }else{
    env <- NULL
    beta_env <- NULL
  }

  # Survey catchability
  Q <- getCatchability(df.tmb, sas)
  Q <- array(Q$Q, dim = c(df.tmb$nage, df.tmb$nsurvey), )
  Q[is.na(Q)] <- 0




  df.fmsy <- list(
    years = years,
    nseason = df.tmb$nseason,
    age = df.tmb$age,
    nage = length(df.tmb$age),
    F0 = matrix(0, length(years), df.tmb$nseason),
    Fsel = Fsel,
    propM = propM,
    propF = propF,
    M = M,
    mat = mat,
    weca = weca,
    west = west,
    sel = Fsel,
    betaSR = betaSR,
    Fbarage = df.tmb$Fbarage,
    nsurvey = df.tmb$nsurvey,
    surveySeason = df.tmb$surveySeason,
    surveyStart = df.tmb$surveyStart,
    surveyEnd = df.tmb$surveyEnd,
    surveySD = 0,
    Q = Q,
    recruitment = recruitment,
    rseason = df.tmb$recseason,
    Fmodel = "sim",
    Ninit = Ninit,
    Rin = NA,
    move = FALSE,
    R0 = R0,
    h = exp(h_val), # Fix this later. need to find it in 'sas'
    logSDR = log(0),
    beta_env = beta_env,
    logSDcatch = log(0),
    env = env,
    alpha = alphaSR,
    b = rep(0, length(years))
  )


  tmp0 <- run.agebased.sms.op(df.fmsy)

  Fmsyin <- seq(0.01 / max(Fselin$selec), 1/ max(Fselin$selec), length.out = 50) # scale to selectivity

  if(stochastic == TRUE){
    df.fmsy$logSDR <- as.numeric(sas$opt$par['logSDrec'])
    message('Running stochastic simulations may take a while')
  }else{
    nruns <- 1
  }

  for (i in 1:length(Fmsyin)) {
    df.fmsy$F0 <- matrix(Fmsyin[i], length(years), df.tmb$nseason)
    for(k in 1:nruns){

    tmp <- run.agebased.sms.op(df.fmsy)

    if (i == 1 & k == 1) {
      df.out <- data.frame(F0 = tmp$Fbar[length(tmp$Fbar)],
                           catch = tmp$Catch[length(years)],
                           SSB = tmp$SSB[length(years)])
    } else {
      df.out <- rbind(
        df.out,
        data.frame(F0 = tmp$Fbar[length(tmp$Fbar)],
                   catch = tmp$Catch[length(years)],
                   SSB = tmp$SSB[length(years)])
      )
    }
  }
}

  if (plotMSY == TRUE) {

    if(stochastic == FALSE){


      print(ggplot(df.out, aes(x = F0, y = catch)) +
      geom_line() +
      theme_classic() +
      scale_x_continuous("Fbar"))
    }else{

      df.mean <- df.out %>% dplyr::group_by(F0) %>%
        dplyr::summarise(catch_m = median(catch),
                  ssb_m = median(SSB))

      # Remove the tail end
      idx.rm <- df.out$F0[which(df.out$SSB < (0.1 * max(df.out$SSB)))]


      print(ggplot(df.out, aes(x = F0, y = catch)) +
              geom_point(alpha = 0.05, color = 'gray') +
              geom_line(data = df.mean, aes(x= F0, y = catch_m), linewidth = 1.3)+
              theme_classic() +
              scale_x_continuous("Fbar"))
    }
  }

if(stochastic == FALSE){
  Fmsy <- list(
    Fmsy = df.out$F0[which.max(df.out$catch)],
    MSY = max(df.out$catch),
    Bmsy = df.out$SSB[which.max(df.out$catch)]
  )
}else{
  df.mean <- df.out %>% dplyr::group_by(F0) %>%
    dplyr::summarise(catch_m = median(catch),
                     ssb_m = median(SSB))

  Fmsy <- list(
    Fmsy = df.mean$F0[which.max(df.mean$catch_m)],
    MSY = max(df.mean$catch_m),
    Bmsy = df.mean$ssb_m[which.max(df.mean$catch_m)]
  )

}

  return(Fmsy)
}
