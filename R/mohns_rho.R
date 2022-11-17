#' Calculate mohns rho.
#'
#' @param df.tmb list of tmb parameters
#' @param peels number of peels
#' @param parms parameters that are estimated
#' @param mps parameters that are mapped
#' @param lower lower boundary of parameters
#' @param upper upper boundary of parameters
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
                      mps = NULL,
                      lwr = list(NA),
                      upr = list(NA),
                      plotfigure = TRUE,
                      limits = NULL,
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


  SSB.base <- getSSB(df.tmb, asses1)
  recruit.base <- getR(df.tmb, asses1)


  F0base <- getFbar(df.tmb, asses1)


  df.save <- data.frame(years = df.tmb$years,
                        SSB = SSB.base$SSB[1:(nrow(SSB.base)-1)],
                        R = recruit.base$R[1:(nrow(SSB.base)-1)],
                        Fbar = F0base$Fbar,
                        peel = 0,
                        convergence = asses1$reps$pdHess)

  df.var <- data.frame(years = df.tmb$years,
                       SSBmin = SSB.base$low[1:(nrow(SSB.base)-1)],
                       SSBmax = SSB.base$high[1:(nrow(SSB.base)-1)],
                       Rmin = recruit.base$low[1:(nrow(SSB.base)-1)],
                       Rmax = recruit.base$high[1:(nrow(SSB.base)-1)],
                       Fmin = F0base$low,
                       Fmax = F0base$high,
                       peel = 0,
                       convergence = 1)


  xx <- data.frame(base = df.save$R,
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

  assess.new <- runAssessment(df.new, parms = parms.new, mps = mps.new)


  SSB.tmb <- getSSB(df.new, assess.new)
  recruit.tmb <- getR(df.new, assess.new)


  F0.tmb <- getFbar(df.new, assess.new)



  tmp <- data.frame(years = df.new$years,
                                 SSB = SSB.tmb$SSB[1:(nrow(SSB.tmb)-1)],
                                 R = recruit.tmb$R[1:(nrow(SSB.tmb)-1)],
                                 Fbar = F0.tmb$Fbar,
                                 peel = i,
                                 convergence = assess.new$reps$pdHess)

  xx[1:(df.new$nyears),i+1] <- tmp$R

  df.save <- rbind(df.save, tmp)

  # Calculate mohns rho
  if(i == 1){
    mohns <- data.frame(
                      SSB = (tmp$SSB[df.new$nyears]-SSB.base$SSB[df.tmb$nyears-i])/SSB.base$SSB[df.tmb$nyears-i],
                      R = (tmp$R[df.new$nyears]-recruit.base$R[df.tmb$nyears-i])/recruit.base$R[df.tmb$nyears-i],
                      F0 = (tmp$Fbar[df.new$nyears]-F0base$Fbar[df.tmb$nyears-i])/F0base$Fbar[df.tmb$nyears-i],
                      row.names = FALSE
                      )
  }else{
    tmp.m <- data.frame(
      SSB = (tmp$SSB[df.new$nyears]-SSB.base$SSB[df.tmb$nyears-i])/SSB.base$SSB[df.tmb$nyears-i],
      R = (tmp$R[df.new$nyears]-recruit.base$R[df.tmb$nyears-i])/recruit.base$R[df.tmb$nyears-i],
      F0 = (tmp$Fbar[df.new$nyears]-F0base$Fbar[df.tmb$nyears-i])/F0base$Fbar[df.tmb$nyears-i],
      row.names = FALSE
    )

    mohns <-  rbind(tmp.m, mohns)
  }

}
  # Plot the figures and calculate Mohns rho

  mohns.tot <- mohns %>% dplyr::summarise(SSB = mean(SSB),
                                   R = mean(R),
                                   F0 = mean(F0))

  df.save <- df.save %>% dplyr::arrange(years)


  df.plot <- df.save #%>% tidyr::pivot_longer(c(SSB,R,Fbar))

  SSBname <- paste('SSB, rho = ', round(mohns.tot$SSB,3))
  Fbarname <- paste('Fbar, rho = ', round(mohns.tot$F0,3))
  Rname <- paste('R, rho = ', round(mohns.tot$R,3))
  #
  #
  if(is.null(limits)){
    limits <- c(min(df.plot$years), max(df.plot$years))
  }

  p1 <- function(){

    x1 <- ggplot2::ggplot(df.plot, aes(x = years, y = SSB, color = factor(peel)))+geom_line()+
      theme_classic()+theme(legend.position = 'none')+
      scale_y_continuous('')+coord_cartesian(xlim = limits)+
      geom_ribbon(data = df.var, aes(ymin = SSBmin, ymax = SSBmax, y = 0), fill = 'red', alpha = .1,
                  linetype = 0)+
      scale_x_continuous('')+
      labs(title=SSBname)

    x2 <- ggplot2::ggplot(df.plot, aes(x = years, y = R, color = factor(peel)))+geom_line()+
      theme_classic()+theme(legend.position = 'none')+
      scale_y_continuous('')+coord_cartesian(xlim = limits)+
      geom_ribbon(data = df.var, aes(ymin = Rmin, ymax = Rmax, y = 0), fill = 'red', alpha = .1,
                  linetype = 0)+
      scale_x_continuous('')+
      labs(title=Rname)

    x3 <- ggplot2::ggplot(df.plot, aes(x = years, y = Fbar, color = factor(peel)))+geom_line()+
      theme_classic()+theme(legend.position = 'none')+
      scale_y_continuous('')+coord_cartesian(xlim = limits)+
      geom_ribbon(data = df.var, aes(ymin = Fmin, ymax = Fmax, y = 0), fill = 'red', alpha = .1,
                  linetype = 0)+
      labs(title=Fbarname)


    p1 <- gridExtra::grid.arrange(x1,x2,x3, nrow =3)

    return(p1)
    }

  if(plotfigure == TRUE){
  print(p1())
  }


  return(list(df.save = df.save,
              p1 = p1(),
              mohns = mohns.tot))
}
