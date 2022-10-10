#' Title
#'
#' @param df.tmb list of tmb parameters
#' @param sas estimated model.
#' @param Fbarage Fbar age
#'
#' @return
#' @export
#'
#' @examples
#' smsplots(df.tmb, reps)
#'
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_classic geom_line geom_ribbon
#' @importFrom ggplot2 alpha scale_y_continuous coord_cartesian theme
#'
 smsPlots <- function(df.tmb, sas, Fbarage = df.tmb$age[df.tmb$age>0]){


  # Plot SSB

  # ages = 0:4
  # nseason = 2
  #

  SSB <- getSSB(df.tmb, sas)

  pssb <- function(){
    ggplot(SSB, aes(x = years, y = SSB))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('SSB')+theme(legend.position = c(0.8,.8))
  }
  # pssb
  # Recruitment


  rec <- getR(df.tmb,sas)

  # take care of crazy recruitment stuff
  lims <- c(min(rec$R)/2, max(rec$R)*2)

  prec <- function(){
    ggplot(rec, aes(x = years, y = R))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('recruitment')+theme(legend.position = 'none')+coord_cartesian(ylim = lims)
  }
  # Fishing mortality

  Fbar <- getFbar(df.tmb, sas, Fbarage = Fbarage)

  pF0 <- function(){
    ggplot(Fbar, aes(x = years, y = Fbar))+geom_line(size = 1.3)+theme_classic()+
    geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('fishing mortality')+theme(legend.position = 'none')
  }

  # And catch
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)
  years <- df.tmb$years
  # Plot SSB, recruitment, catch and fishing mortality

  Catch <- getCatch(df.tmb,sas)

  pCatch <- function(){
    ggplot(Catch, aes(x = years, y = Catch))+geom_line(size = 1.4)+
    theme_classic()+geom_ribbon(aes(ymin = minSE, ymax = maxSE), fill = alpha('red', 0.2), linetype = 0)+
    scale_y_continuous('Catch')+theme(legend.position = c(0.8,.8))
  }

  ls <- gridExtra::grid.arrange(pssb(),pCatch(),prec(),  pF0(), ncol = 2)

  ls

  return(ls)
  }
