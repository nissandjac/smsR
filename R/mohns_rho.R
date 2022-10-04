#' Calculate mohns rho.
#'
#' @param df.tmb list of tmb parameters
#' @param peels number of peels
#' @param parms parameters that are estimated
#' @param mps parameters that are mapped
#' @param lower lower boundary of parameters
#' @param upper upper boundary of parameters
#' @param Fbarage ages to calculate Fbar
#' @param plotfigure TRUE FALSE of whether to plot figure
#' @param limits limits on plot
#' @param dll which dll to use. Uses 'sms' as standard
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_classic geom_line geom_ribbon
#' @importFrom ggplot2 alpha scale_y_continuous coord_cartesian theme facet_wrap

mohns_rho <- function(df.tmb,
                      peels = 5,
                      parms ,
                      mps,
                      lwr = list(NA),
                      upr = list(NA),
                      Fbarage = c(1,2),
                      plotfigure = TRUE,
                      limits = NA,
                      dll = 'smsR'
                      ){

  # Remove the dplyr summarise warning
  options(dplyr.summarise.inform = FALSE)

  # Run first the model as is
  # obj <-TMB::MakeADFun(df.tmb,parms,DLL=dll, map = mps, silent = TRUE)
  #
  # x <- obj$report()
  # system.time(opt<-stats::nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
  #                         control = list(iter.max = 1e6,
  #                                        eval.max = 1e6))) #
  # system.time(rep<-TMB::sdreport(obj))
  # rep

  # Save the results to a data frame
  asses1 <- runAssessment(df.tmb, parms = parms, mps = mps)


  SSB.base <- asses1$reps$value[names(asses1$reps$value) == 'SSB']
  recruit.base <- asses1$reps$value[names(asses1$reps$value) == 'Rsave']


  F0base <- data.frame(F0 = asses1$reps$value[names(asses1$reps$value) == 'F0'])
  F0base$age <- df.tmb$age
  F0base$season <- rep(1:df.tmb$nseason, each = df.tmb$nage*df.tmb$nyears)
  F0base$year <- rep(years, each = length(ages))


  Fbarbase <- F0base[F0base$age %in% Fbarage,] %>% dplyr::group_by(year, age) %>%
    dplyr::summarise(Fbar0 = sum(F0)) %>% dplyr::group_by(year) %>%
    dplyr::summarise(Fmean = mean(Fbar0))



  df.save <- data.frame(years = df.tmb$years,
                        SSB = SSB.base[1:(length(SSB.base)-1)],
                        R = recruit.base[1:(length(SSB.base)-1)],
                        Fbar = Fbarbase$Fmean,
                        peel = 0,
                        convergence = asses1$reps$pdHess)



  xx <- data.frame(base = df.save$SSB,
                   '-1' = NA,
                   '-2' = NA,
                   '-3' = NA,
                   '-4' = NA,
                   '-5' = NA,
                   row.names = df.tmb$years
                   )

  #colnames(xx) <- names(shake)

  for(i in 1:peels){


  df.new <- df.tmb

  if(i == 1){

    rmidx <- length(df.new$years)
    }else{
    rmidx <- ((length(df.new$years)-(i-1)):length(df.new$years))
  }

  df.new$years <- df.new$years[-rmidx]

  df.new$nyears <- length(df.new$years)
  df.new$effort <- df.new$effort[1:df.new$nyears,, drop = F]
  df.new$bidx <- df.new$bidx[1:df.new$nyears]
  df.new$Surveyobs <- df.new$Surveyobs[,1:df.new$nyears,, drop = F]
  df.new$Catchobs <- df.new$Catchobs[,1:df.new$nyears,, drop = F]
  df.new$nocatch <- df.new$nocatch[1:df.new$nyears,, drop = F]


  parms.new <- parms
  parms.new$logRin <- parms.new$logRin[1:df.new$nyears]
  parms.new$Fyear <- parms.new$logRin[1:df.new$nyears]

  if(max(df.new$bidx) > 0){
    if(max(df.new$years) < max(df.new$years[df.new$bidx])){
      stop('peeling below the latest selectivity year')
    }
  }

  mps.new <- mps

  if('Fyear' %in% names(mps.new)){
    mps.new$Fyear <- factor(rep(NA, df.new$nyears))
  }

  if('logRin' %in% names(mps.new)){
    mps.new$logRin <- factor(rep(NA, df.new$nyears))
  }

  assess.new <- runAssessment(df.new, parms = parms.new, mps = mps.new, lwr = lwr, upr = upr)

  SSB.tmb <- assess.new$reps$value[names(assess.new$reps$value) == 'SSB']
  recruit.tmb <- assess.new$reps$value[names(assess.new$reps$value) == 'Rsave']


  F0 <- data.frame(F0 = assess.new$reps$value[names(assess.new$reps$value) == 'F0'])
  F0$age <- df.new$age
  F0$season <- rep(1:df.new$nseason, each = df.new$nage*df.new$nyears)
  F0$year <- rep(df.new$years, each = length(df.new$age))


  Fbar <- F0[F0$age %in% Fbarage,] %>% dplyr::group_by(year, age) %>%
    dplyr::summarise(Fbar0 = sum(F0)) %>% dplyr::group_by(year) %>%
    dplyr::summarise(Fmean = mean(Fbar0))


  tmp <- data.frame(years = df.new$years,
                                 SSB = SSB.tmb[1:(length(SSB.tmb)-1)],
                                 R = recruit.tmb[1:(length(SSB.tmb)-1)],
                                 Fbar = Fbar$Fmean,
                                 peel = i,
                                 convergence = assess.new$reps$pdHess)

  xx[1:(df.new$nyears),i+1] <- tmp$SSB

  df.save <- rbind(df.save, tmp)

  # Calculate mohns rho
  if(i == 1){
    mohns <- data.frame(
                      SSB = (SSB.tmb[df.new$nyears]-SSB.base[df.tmb$nyears-i])/SSB.base[df.tmb$nyears-i],
                      R = (tmp$R[df.new$nyears]-recruit.base[df.tmb$nyears-i])/recruit.base[df.tmb$nyears-i],
                      F0 = (tmp$Fbar[df.new$nyears]-Fbarbase$Fmean[df.tmb$nyears-i])/Fbarbase$Fmean[df.tmb$nyears-i],
                      row.names = FALSE
                      )
  }else{
    tmp.m <- data.frame(
      SSB = (tmp$SSB[df.new$nyears]-SSB.base[df.tmb$nyears-i])/SSB.base[df.tmb$nyears-i],
      R = (tmp$R[df.new$nyears]-recruit.base[df.tmb$nyears-i])/recruit.base[df.tmb$nyears-i],
      F0 = (tmp$Fbar[df.new$nyears]-Fbarbase$Fmean[df.tmb$nyears-i])/Fbarbase$Fmean[df.tmb$nyears-i],
      row.names = FALSE
    )

    mohns <-  rbind(tmp.m, mohns)
  }

}
  # Plot the figures and calculate Mohns rho

  mohns.tot <- mohns %>% dplyr::summarise(SSB = mean(SSB),
                                   R = mean(R),
                                   F0 = mean(F0))



  df.plot <- df.save %>% tidyr::pivot_longer(2:4)

  df.plot$name[df.plot$name == 'SSB'] <- paste('SSB, rho = ', round(mohns.tot$SSB,3))
  df.plot$name[df.plot$name == 'Fbar'] <- paste('Fbar, rho = ', round(mohns.tot$F0,3))
  df.plot$name[df.plot$name == 'R'] <- paste('R, rho = ', round(mohns.tot$R,3))


  if(is.na(limits)){
    limits <- c(min(df.plot$years), max(df.plot$years))
  }

  p1 <- function(){
    ggplot2::ggplot(df.plot, aes(x = years, y = value, color = factor(peel)))+geom_line()+
    ggplot2::facet_wrap(~name, scales = 'free_y', nrow = 3)+theme_classic()+theme(legend.position = 'none')+
    scale_y_continuous('')+coord_cartesian(xlim = limits)
  }

  if(plotfigure == TRUE){
  print(p1())
  }


  return(list(df.save = df.save,
              p1 = p1(),
              mohns = mohns.tot))
}
