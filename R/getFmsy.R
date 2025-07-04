#' Calculate Fmsy from a sms model
#'
#' @param df.tmb sms data frame
#' @param sas Fitted sms model
#' @param nyears number of years to simulate
#' @param recruitment recruitment type to calculate Fmsy. Currently only hockey works.
#' @param plotMSY TRUE OR FALSE to plot the result
#' @param stochastic stochastic fmsy?
#'
#' @return returns the fishing mortality that leads to Fmsy
#' @export
#'
getFmsy <- function(sas,
                    nyears = 50,
                    recruitment = "hockey",
                    plotMSY = FALSE,
                    stochastic = FALSE) {
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

  mat <- array(df.tmb$Mat[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  weca <- array(df.tmb$weca[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  west <- array(df.tmb$west[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  propM <- array(df.tmb$propM[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))
  propF <- array(df.tmb$propF[, df.tmb$nyears, ], dim = c(df.tmb$nage, nyears, df.tmb$nseason))

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
    df.fmsy$F0 <- matrix(Fmsyin[i], length(years), df.tmb$nseason)

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
