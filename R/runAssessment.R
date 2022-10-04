#' Run the sms model
#'
#' @param df.tmb list of tmb parameters
#' @param lower lower bounds
#' @param upper upper bounds
#' @param parms optional parameter input
#' @param mps optional mapping input
#'
#' @return
#' a list of derived variables
#' @export
#'
#' @examples
#' x <- runAssessment(df.tmb)
runAssessment <- function(df.tmb,
                          lower = -Inf,
                          upper = Inf,
                          parms = NA,
                          mps = NA){





# Try initial parameters in weird spots
SDCatch <- array(rep(.6,length(unique(df.tmb$catchCV))), dim = c(length(unique(df.tmb$catchCV)),1))
SDsurvey <- rep(.1, length(unique(df.tmb$Qidx_CV[df.tmb$Qidx_CV > -1])))


if(is.list(parms)){
  # Check lengths here
}else{

  parms <- list(logRin = rep(18, df.tmb$nyears),
              logNinit = rep(15, df.tmb$nage-1),
              Fyear = rep(1, df.tmb$nyears), # Mapped out
              Fseason = matrix(0.5, nrow = 1, ncol = length(unique(df.tmb$bidx))),
              logFage = matrix(1, nrow = df.tmb$Fmaxage+1, ncol = length(unique(df.tmb$bidx))),
              SDsurvey = SDsurvey,
              SDcatch = SDCatch,
              logQ = rep(log(1), 5),#length(df.tmb$surveyCV)
              pin = 1,
              logalpha = 2,
              logbeta = df.tmb$betaSR,
              logSDrec = log(0.5))
}
# Not estimated parameters
mps <- list(
  # 'pin' = factor(parms$pin*NA),
  'Fyear' = parms$Fyear,
  'SDcatch' = factor(parms$SDcatch*NA),
  'logbeta' = factor(parms$logbeta*NA)
)



obj <-TMB::MakeADFun(df.tmb,parms,DLL="sms", map = mps)
x <- obj$report()

system.time(opt<-stats::nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
                        control = list(iter.max = 1e6,
                                       eval.max = 1e6))) #


system.time(rep<-TMB::sdreport(obj))


return(list(x = x,
            opt = opt,
            rep = rep))
}
