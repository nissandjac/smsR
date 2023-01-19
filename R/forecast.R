
#' Forecast an sms model
#'
#' @param df.tmb sms dataframe
#' @param sas fitted sms model
#' @param recruitment hockey, long_mean or short_mean
#' @param HCR Fmsy or Bescape
#' @param avg_R years to average R
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
#                     avg_years = rep(3, 4),
                     avg_R = NULL,
                     Bpa = NULL,
                     Fcap = NULL
                     ){





  # Extract N in the beginning of the coming year

  # Get the estimated survey
  reps <- sas$reps
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)

  N_temp <- as.numeric(sdrep[rep.values == 'term_logN_next',1])

  Rold <- getR(df.tmb, sas)[-df.tmb$nyears+1,]

  if(is.null(avg_R)){
    avg_R <- df.tmb$years
  }

  if(recruitment == 'mean'){
    N_temp[1] <- log(exp(mean(log(Rold$R[Rold$years %in% avg_R]), na.rm = TRUE)))
  }

  N_temp <- exp(N_temp)


  if(df.tmb$nseason == 1){
    Fsel <- getF(df.tmb, sas)
    Fsel <- Fsel$F0[Fsel$years == max(df.tmb$years)]
    Fsel[df.tmb$age < df.tmb$Fminage & df.tmb$age > df.tmb$Fmaxage] <- 0

    Fsel <- (Fsel/mean(Fsel[(df.tmb$Fbarage[1]:df.tmb$Fbarage[2]+1)])) # Scale selectivity to Fbar
  }else{

    N_current <- matrix(0 , nrow = df.tmb$nage, ncol = df.tmb$nseason)
    N_current[2:df.tmb$nage,1] <- N_temp[2:df.tmb$nage]
    N_current[1, df.tmb$recseason] <- N_temp[1]

    Fsel <- getF(df.tmb, sas)
    Fsel <- matrix(Fsel$F0[Fsel$years == max(df.tmb$years)], nrow = df.tmb$nage, ncol = df.tmb$nseason)
    Fsel[df.tmb$age < df.tmb$Fminage & df.tmb$age > df.tmb$Fmaxage] <- 0


    Fsel <- (Fsel/mean(rowSums(Fsel)[(df.tmb$Fbarage[1]:df.tmb$Fbarage[2])+1])) # Scale selectivity to Fbar
  }
  # Now do a forecast with fishing mortality

  if(HCR == 'Fmsy'){ # Calc Fmsy from assessment model

    Fmsy <- getFmsy(df.tmb, sas, recruitment = 'hockey')
    F0 <- matrix(Fsel*Fmsy$Fmsy, nrow = df.tmb$nage, ncol = df.tmb$nseason)
    fc <- forecast.sms(df.tmb, N_current, F0)

  }

  if(HCR == 'Bescape'){

    if(is.null(Fcap)){warning('assuming Fcap = 2yr-1')}

    # Find the fishing mortality leading to Bescape
    F0 <- calcBescape(Bpa, df.tmb, N_current, Fsel, Fcap)*Fsel
    # Do forecast
    fc <- forecast.sms(df.tmb , N_current , F0)
    Fmsy <- mean(rowSums(F0)[(df.tmb$Fbarage[1]:df.tmb$Fbarage[2])+1])


    Fmsy <- list(Fmsy = Fmsy, MSY = NA)


  }







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
forecast.sms <- function(df.tmb , N_current, F0 ){

  N_future <- matrix(NA, df.tmb$nage) # For the following year SSB
  C_new <- matrix(0, df.tmb$nage, df.tmb$nseason)

  #
  M <- df.tmb$M[,df.tmb$nyears+1,]
  weca <- df.tmb$weca[,df.tmb$nyears+1,]
  west <- df.tmb$west[,df.tmb$nyears+1,]
  mat <- df.tmb$Mat[,df.tmb$nyears+1,]



  Z <- F0 + M


  for (qrts in 1:df.tmb$nseason) {
    for (i in 1:df.tmb$nage){

      if(Z[i,qrts]>0){
        C_new[i,qrts]=(F0[i,qrts]/Z[i,qrts])*N_current[i,qrts]*weca[i,qrts]*(1-exp(-(Z[i,qrts])))
      }

      if(qrts < df.tmb$nseason){

        if(i == 1 & qrts < df.tmb$recseason){
         # Dont do anything
        }else{
         N_current[i, qrts+1] <- N_current[i,qrts]*exp(-Z[i, qrts])
        }
      }else{

        N_future[1] <- 0
        if(i < df.tmb$nage){
          N_future[i+1] <- N_current[i,qrts]*exp(-Z[i,qrts])
        }else{
          N_future[df.tmb$nage] <- N_future[df.tmb$nage]+N_current[i,qrts]*exp(-Z[i,qrts])
        }
      }
        if(mat[1,1]>0){
          warning('new recruits not included in SSB for forecasted year')
        }
         SSB_next <- sum(N_future*weca[,1]*mat[,1])
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
calcFTAC <- function(TAC , df.tmb){


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




  data.in <- list(TAC = TAC,
                  Ncurrent = Ncurrent,
                  df.tmb = df.tmb)

  parms.in <- list(0.5)

  Fnew <- optim(parms.in, lower = 0.0001, upper = 2, fn = optTAC, data= data.in, method = 'L-BFGS-B')



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
calcBescape <- function(Bpa ,
                        df.tmb,
                        N_current,
                        Fsel,
                        Fcap = 2){

  optBpa <- function(data, par ){
    df.tmb <- data$df.tmb
    Fsel <- data$Fsel
    Fcalc <- as.numeric(par[1])
    Bpa <- data$Bpa

    F0 <- Fsel*Fcalc
    ls <- forecast.sms(df.tmb, N_current,F0)

    ans <- (Bpa - ls$SSB)^2
  }


  data.in <- list(Bpa = Bpa,
                  N_current = N_current,
                  Fsel = Fsel,
                  df.tmb = df.tmb)

  parms.in <- list(0.5)

  Fnew <- optim(parms.in, lower = 0.0001, upper = Fcap, fn = optBpa, data= data.in, method = 'L-BFGS-B')

  return(Fnew$par)
}




#' Produce an ICES forecast table
#'
#' @param df.tmb sms model input data
#' @param sas fitted sms model
#' @param TACold last years tac
#' @param HCR Harvest control rule to use - options Bescape or Fmsy
#'
#' @return
#' @export
#'
#' @examples
createForecastTable <- function(df.tmb, sas, TACold, HCR = 'Bescape', avg_R = 'mean'){

  # Get numbers at age
  # Get the estimated survey
  reps <- sas$reps
  sdrep <- summary(reps)
  rep.values<-rownames(sdrep)

  N_temp <- as.numeric(sdrep[rep.values == 'term_logN_next',1])

  Rold <- getR(df.tmb, sas)[-df.tmb$nyears+1,]

  if(is.null(avg_R)){
    avg_R <- df.tmb$years
  }

  if(recruitment == 'mean'){
    N_temp[1] <- log(exp(mean(log(Rold$R[Rold$years %in% avg_R]), na.rm = TRUE)))
  }

  N_temp <- exp(N_temp)

  if(df.tmb$nseason == 1){
    Fsel <- getF(df.tmb, sas)
    Fsel <- Fsel$Fsel[Fsel$years == max(df.tmb$years)]
    Fsel[df.tmb$age < df.tmb$Fminage & df.tmb$age > df.tmb$Fmaxage] <- 0

    Fsel <- (Fsel/mean(Fsel[(df.tmb$Fbarage[1]:df.tmb$Fbarage[2])+1])) # Scale selectivity to Fbar
  }else{

    N_current <- matrix(0 , nrow = df.tmb$nage, ncol = df.tmb$nseason)
    N_current[2:df.tmb$nage,1] <- N_temp[2:df.tmb$nage]
    N_current[1, df.tmb$recseason] <- N_temp[1]

    Fsel <- getF(df.tmb, sas)
    Fsel <- matrix(Fsel$F0[Fsel$years == max(df.tmb$years)], nrow = df.tmb$nage, ncol = df.tmb$nseason)
    Fsel[df.tmb$age < df.tmb$Fminage & df.tmb$age > df.tmb$Fmaxage] <- 0


    Fsel <- (Fsel/mean(rowSums(Fsel)[(df.tmb$Fbarage[1]:df.tmb$Fbarage[2])+1])) # Scale selectivity to Fbar
  }


  F0 <- matrix(0, df.tmb$nage, df.tmb$nseason)
  f0 <- forecast.sms(df.tmb , N_current , F0)

  # Forecast with last years fishing
  Flast <- getF(df.tmb, sas) %>% dplyr::filter(years == max(df.tmb$years))
  F0 <- matrix(Flast$F0, nrow = df.tmb$nage, ncol = df.tmb$nseason)

  forlast <- forecast.sms(df.tmb, N_current, F0)

  # Bpa = Blim
  Flim <- calcBescape(Bpa = df.tmb$betaSR, df.tmb, N_current, Fsel, Fcap)*Fsel
  # Do forecast
  fc <- forecast.sms(df.tmb , N_current , F0)
  Fmsy <- mean(rowSums(F0)[df.tmb$Fbarage[1]:df.tmb$Fbarage[2]])




}





