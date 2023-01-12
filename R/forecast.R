
#' Forecast an sms model
#'
#' @param df.tmb sms
#' @param sas
#' @param recruitment
#' @param HCR
#' @param avg_years years to average weca, M, west, and mat
#' @param method
#'
#' @return
#' @export
#'
#' @examples
getTAC <- function(df.tmb,
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
    Fsel <- x$Fsel[x$years == max(df.tmb$years)]
    Fsel[Fsel < df.tmb$Fminage & Fsel > df.tmb$Fmaxage] <- 0

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

  fc <- forecast.sms(df.tmb, F0, M , N_current, weca, west)




ls.out <- list(TAC = sum(C_new),
               SSB = SSB_next,
               N_future = N_future,
               Fmsy = Fmsy$Fmsy,
               MSY = Fmsy$MSY)


    }




#' Simple forecast of one year.
#'
#' @param df.tmb
#' @param F0
#' @param M
#' @param N_current
#' @param weca
#' @param west
#'
#' @return
#' @export
#'
#' @examples
forecast.sms <- function(df.tmb, F0, M, N_current, weca, west){

  N_new <- matrix(NA, df.tmb$nage, df.tmb$nseason)
  N_future <- matrix(NA, df.tmb$nage) # For the following year SSB
  C_new <- matrix(NA, df.tmb$nage, df.tmb$nseason)

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
            SSB = sum(SSB_next),
            )
       )
}
