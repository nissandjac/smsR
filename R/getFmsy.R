#' Calculate Fmsy from a sms model
#'
#' @param sas Fitted sms model
#' @param nyears number of years to simulate
#' @param recruitment recruitment type to calculate Fmsy. Currently only hockey works.
#' @param plotMSY TRUE OR FALSE to plot the result
#' @param stochastic stochastic fmsy?
#' @param nruns number of replicates in stochastic fmsy
#'
#' @return returns the fishing mortality that leads to Fmsy
#' @export
#'
getFmsy <- function(sas,
                    nyears = 50,
                    recruitment = "hockey",
                    plotMSY = FALSE,
                    stochastic = FALSE,
                    nruns = 1000) {
  #  Number
  df.tmb <- sas$dat
  years <- 1:nyears
  nyears <- length(years)

  if(sas$dat$randomM == 1){

  M_estimated <- getM_var(sas$dat, sas) %>% filter(years %in% max(sas$dat$years))


  M <-   array(M_estimated$M_var, dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  }else{
  M <- array(df.tmb$M[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  }

  mat_last <- df.tmb$Mat[, df.tmb$nyears, , drop = FALSE]   # nage × 1 × nseason
  mat <- mat_last[, rep(1, nyears), , drop = FALSE]


  weca_last <- df.tmb$weca[, df.tmb$nyears, , drop = FALSE]   # nage × 1 × nseason
  weca <- weca_last[, rep(1, nyears), , drop = FALSE]

  west_last <- df.tmb$west[, df.tmb$nyears, , drop = FALSE]   # nage × 1 × nseason
  west <- west_last[, rep(1, nyears), , drop = FALSE]

  propM_last <- df.tmb$propM[, df.tmb$nyears, , drop = FALSE]   # nage × 1 × nseason
  propM <- propM_last[, rep(1, nyears), , drop = FALSE]

  propF_last <- df.tmb$propF[, df.tmb$nyears, , drop = FALSE]   # nage × 1 × nseason
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
    alphaSR <- exp(parms.est$value[parms.est$parameter == "logalpha"])
    betaSR <- df.tmb$betaSR
    R0 <- alphaSR * df.tmb$betaSR
  } else {
    alphaSR <- NA
    betaSR <- NA
    R0 <- exp(parms.est$value[parms.est$parameter == "logR0"])
    h_val <- attr(sas$obj$env$parameters$logh, "shape")

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
    Ninit = NULL,
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
