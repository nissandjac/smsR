
#' Forecast an sms model
#'
#' @param df.tmb sms dataframe
#' @param sas fitted sms model
#' @param recruitment hockey, long_mean or short_mean
#' @param HCR Fmsy or Bescape
#' @param avg_years years to average weca, M, west, and mat
#' @param method
#'
#' @return
#' @export
#'
#' @examples
calcTAC <- function(df.tmb,
                     sas,
                     recruitment = 'hockey',
                     HCR = 'Fmsy',
                     avg_years = rep(3, 4)
                     ){




  if(length(avg_years) == 1){
    avg_years <- rep(avg_years, 4)
  }
  # Extract N in the beginning of the coming year

  # Get the estimated survey
  reps <- sas$reps
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)

  N_current <- as.numeric(sdrep[rep.values == 'term_logN_next',1])

  Rold <- getR(df.tmb, sas)


  if(recruitment == 'long_mean'){
    N_current[1] <- log(mean(Rold$R, na.rm = TRUE))
  }
  if(recruitment == 'short_mean'){
    N_current[1] <- log(mean(Rold$R[Rold$years %in% df.tmb$years[(df.tmb$nyears-5):df.tmb$nyears]]))
  }

  N_current <- exp(N_current)


  if(df.tmb$nseason == 1){
    Fsel <- getSelex(df.tmb, sas)
    Fsel <- Fsel$Fsel[Fsel$years == max(df.tmb$years)]
    Fsel[df.tmb$age < df.tmb$Fminage & df.tmb$age > df.tmb$Fmaxage] <- 0

    Fsel <- (Fsel/mean(Fsel[df.tmb$Fbarage[1]:df.tmb$Fbarage[2]])) # Scale selectivity to Fbar
  }else{
    stop('add this feature for > 1 season')
  }
  # Now do a forecast with fishing mortality

  if(HCR == 'Fmsy'){ # Calc Fmsy from assessment model

    Fmsy <- getFmsy(df.tmb, sas, recruitment = 'hockey')
    F0 <- matrix(Fsel*Fmsy$Fmsy, nrow = df.tmb$nage, ncol = df.tmb$nseason)

    }else{

      stop('Fmsy only option right now, get coding for Bescape and others')
  }

  M <- matrix(rowMeans(df.tmb$M[,(df.tmb$nyears-avg_years[1]+1):df.tmb$nyears,]), nrow = df.tmb$nage, ncol = df.tmb$nseason)
  mat <- matrix(rowMeans(df.tmb$Mat[,(df.tmb$nyears-avg_years[2]+1):df.tmb$nyears,]), nrow = df.tmb$nage, ncol = df.tmb$nseason)
  weca <- matrix(rowMeans(df.tmb$weca[,(df.tmb$nyears-avg_years[3]+1):df.tmb$nyears,]), nrow = df.tmb$nage, ncol = df.tmb$nseason)
  west <- matrix(rowMeans(df.tmb$weca[,(df.tmb$nyears-avg_years[4]+1):df.tmb$nyears,]), nrow = df.tmb$nage, ncol = df.tmb$nseason)


  bios <- list(M = M,
               mat = mat,
               weca = weca,
               west = west)


  fc <- forecast.sms(df.tmb, N_current, F0, bios)




ls.out <- list(TAC = fc$Catch,
               SSB = fc$SSB,
               Fmsy = Fmsy$Fmsy,
               Fnext = F0,
               MSY = Fmsy$MSY)


    }




#' Simple one year forecast
#'
#' @param df.tmb
#' @param N_current
#' @param F0
#' @param bios
#'
#' @return
#' @export
#'
#' @examples
forecast.sms <- function(df.tmb ,N_current, F0, bios){

  N_new <- matrix(NA, df.tmb$nage, df.tmb$nseason)
  N_future <- matrix(NA, df.tmb$nage) # For the following year SSB
  C_new <- matrix(NA, df.tmb$nage, df.tmb$nseason)

  #
  M <- bios$M
  weca <- bios$weca
  west <- bios$west
  mat <- bios$mat



  Z <- F0 + M
  N_new[,1] <- N_current


  for (qrts in 1:df.tmb$nseason) {
    for (i in 1:df.tmb$nage){

      C_new[i,qrts]=(F0[i,qrts]/Z[i,qrts])*N_new[i,qrts]*weca[i,qrts]

      if(i < df.tmb$nseason){
        N_new[i, qrts+1] <- N_new*exp(Z[i, qrts])
      }else{
        N_future[1] <- 0
        if(i < df.tmb$nage){
          N_future[i+1] <- N_new[i]*exp(Z[i,qrts])
        }else{
          N_future[i] <- N_future[i,qrts]+N_new[i,qrts]*exp(-Z[i,qrts])
        }
      }
         SSB_next <- sum(N_future*weca*mat)
    }
  }
return(list(Catch = sum(C_new),
            SSB = sum(SSB_next)
            )


       )
}


#' Calculate the fishing mortality required to get the assigned TAC
#'
#' @param TAC TAC were are trying to achieve
#' @param df.tmb SMS input list
#' @param bios list of bioparams
#' @param N_current numbers at age at the beginning of the year
#' @param findTAC Look for TAC or SSB

#'
#' @return
#' @export
#'
#' @examples
calcFTAC <- function(TAC , df.tmb, bios , Ncurrent, findTAC = TRUE){

  if(findTAC == TRUE){
  data.in <- list(bios = bios,
                  TAC = TAC,
                  Ncurrent = Ncurrent,
                  df.tmb = df.tmb)

  parms.in <- list(0.5)

  Fnew <- optim(parms.in, lower = 0.0001, upper = 2, fn = optTAC, data= data.in, method = 'L-BFGS-B')




  }

}


#' Optimizer for finding the fishing mortality that leads to TAC
#'
#' @param data
#' @param par
#'
#' @return
#' @export
#'
#' @examples
optTAC <- function(data, par){

  df.tmb <- data$df.tmb
  Fsel <- data$bios$Fsel
  Fcalc <- as.numeric(par[1])
  TAC <- data$TAC

  F0 <- Fsel*Fcalc


  ls <- forecast.sms(df.tmb, N_current,F0, bios)

  ans <- (log(TAC) - log(ls$Catch))^2

return(ans)
}


