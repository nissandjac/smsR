#' Run the sms model
#'
#' @param df.tmb list of tmb parameters
#' @param lower lower bounds
#' @param upper upper bounds
#' @param parms optional parameter input
#' @param mps optional mapping input
#' @param silent Supress TMB output
#'
#' @return
#' a list of derived variables
#' @export
#'
#' @examples
#' sas <- runAssessment(df.tmb, parms) # Run an assessment
runAssessment <- function(df.tmb,
                          parms,
                          lwr = list(NA),
                          upr = list(NA),
                          mps = NULL,
                          silent = TRUE){


if(is.null(mps)){
  # map the parameters that are not estimated in this model
  mps <- getMPS(df.tmb, parms)

}

if(is.null(df.tmb$betaSR)){
  df.tmb$betaSR <- NA
}

if(length(df.tmb$Qidx) > length(parms$logQ)){
  stop('wrong length of survey ages - redo input parameters')
}

  if(df.tmb$randomF == 1){
    rlist <- c('logFyear','logRin')
  }else{
    rlist <- c()
  }

  rlist <- c()


obj <-TMB::MakeADFun(df.tmb,parms,DLL="smsR", map = mps, silent = silent, random = rlist)
x <- obj$report()

# Set boundaries
lower <- obj$par-Inf
upper <- obj$par+Inf

# Permanent bounds for realism
lower[names(lower) == 'logFyear' ] <- log(0.001)
lower[names(lower) == 'Fseason'] <- 0.0001
lower[names(lower) == 'SDsurvey'] <- 0.0001
lower[names(lower) == 'logSDrec'] <- log(0.1)
lower[names(lower) == 'SDcatch'] <- 0.01
lower[names(lower) == 'creep'] <-  -0.1

upper[names(upper) == 'SDsurvey'] <- 2
upper[names(upper) == 'logSDrec'] <- log(4)
upper[names(upper) == 'logFyear'] <- log(10)
upper[names(upper) == 'SDcatch'] <- sqrt(2)
upper[names(upper) == 'creep'] <- 0.2

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





system.time(opt<-stats::nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
                        control = list(iter.max = 1e6,
                                       eval.max = 1e6))
            ) #


system.time(reps<-TMB::sdreport(obj))


# Add congrats
if(reps$pdHess == TRUE & max(abs(reps$gradient.fixed)) < 0.1){
  message('Congrats, your model converged')
}else{
  message('Check model convergence')
  message(paste('hessian =', reps$pdHess))
  message(max(abs(reps$gradient.fixed)))
}


# Throw a warning if something is on the boundary
if(any(reps$par.fixed == lower)){

  parms.boundary <- reps$par.fixed[reps$par.fixed == lower]

  for(i in 1:length(parms.boundary)){
  warning(paste(names(parms.boundary[i]),'estimated on the lower boundary'))
  }
}

if(any(reps$par.fixed == upper)){

  parms.boundary <- reps$par.fixed[reps$par.fixed == upper]

  for(i in 1:length(parms.boundary)){
    warning(paste(names(parms.boundary[i]),'estimated on the upper boundary'))
  }
}


return(structure(list(x = x,
            opt = opt,
            reps = reps,
            obj = obj), class="sms"))
}
