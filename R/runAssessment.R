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
runAssessment <- function(df.tmb,
                          lwr = list(NA),
                          upr = list(NA),
                          parms,
                          mps){




obj <-TMB::MakeADFun(df.tmb,parms,DLL="smsR", map = mps)
x <- obj$report()

# Set boundaries
lower <- obj$par-Inf
upper <- obj$par+Inf


# Add custom boundaries to parameters
for(i in 1:length(lwr)){

  if(is.na(lwr[[1]][1]) == 0){

  idx <- (names(lower) %in% names(lwr)[i])

  if(sum(idx) != length(lwr[[i]])){

    if(length(lwr[[i]]) != 1){
      stop('Error in custom boundaries')
    }
      }

  lower[idx] <- lwr[[i]]
  }
}
#
for(i in 1:length(upr)){
  if(is.na(upr[[1]][1]) == 0){
    idx2 <- (names(upper) %in% names(upr)[i])
    if(sum(idx2) != length(upr[[i]])){
      if(length(lwr[[i]]) != 1){
      stop('Error in custom boundaries')
      }
    }

    upper[idx2] <- upr[[i]]
  }
}

# Permanent bounds for realism
lower[names(lower) == 'Fyear' ] <- 0.001
lower[names(lower) == 'Fseason'] <- 0.0001
upper[names(upper) == 'SDsurvey'] <- 3
lower[names(lower) == 'logSDrec'] <- log(0.1)
upper[names(upper) == 'logSDrec'] <- log(2)




system.time(opt<-stats::nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
                        control = list(iter.max = 1e6,
                                       eval.max = 1e6))) #


system.time(reps<-TMB::sdreport(obj))


return(list(x = x,
            opt = opt,
            reps = reps))
}
