#' Run a full age based model as operating model for SMS
#'
#' @param df list of parameters and life history values
#'
#' @return
#' @export
#'
#' @examples
#' run.agebased.true.catch(df) # runs the model
#'
run.agebased.sms.op <- function(df){



  nseason <- df$nseason

  df$tEnd <- length(df$years)*nseason
  nyear <-df$tEnd/df$nseason #
  year <- df$years
  tEnd <- nyear*nseason

  # Catchability

  # Maturity
  mat <- df$mat # Fecundity
  mat[is.na(mat)] <- 0


  # Mortality
  M <- df$M
  M[is.na(M)] <- 0 # Remove NA's

  # Weight
  weca <- df$weca
  weca[is.na(weca)] <- 0

  west <- df$west
  west[is.na(west)] <- 0

  # Fishing mortality
  F0 <- df$F0
  F0[is.na(F0)] <- 0
  # Age
  nage <- df$nage
  age <- df$age
  maxage <- max(df$age)


 #Mage <- c(0,cumsum(M[1:(nage-1),1,3]))


  R0 <- df$R0
  SDR <- df$SDR

  # Calculate N0 based on R0
  mage <- max(df$age) # Max age
  agetmp <- 0:(mage*3)
  nagetmp <- length(agetmp)

  M0 <- rep(0, mage*3+1)
  M0[1:nage]<- M[,1,1]*df$nseason
  M0[nage:length(M0)] <- M[nage,1,1]*df$nseason

  N0tmp <- rep(NA,nagetmp)

  N0tmp[1:(nagetmp-1)] = R0*exp(-agetmp[1:(nagetmp-1)]*M0[1:(nagetmp-1)])
  N0tmp[nagetmp] =  R0*exp(-M0[nagetmp]*agetmp[nagetmp])/(1-exp(-M0[nagetmp]))



  N0 <- matrix(NA,nage)
  N0[1:(nage-1)] <- N0tmp[1:(nage-1)]
  N0[nage] <- sum(N0tmp[nage:nagetmp])

  SSB_0 <- sum(N0*west[,1,1]*mat[,1,1])


  # Recruitment per country
  R_0 <- R0
  nspace <- 1

  # Initialize for saving
  year_1 <- c(year,max(year)+1)

  SSB <- matrix(NA,nyear, nspace,
                dimnames = list(year = df$years,
                                space = 1:nspace))
  SSB.all <- array(NA, dim = c(nyear, nspace,nseason),
                   dimnames = list(year = year, space = 1:nspace, season = 1:nseason))
  SSB.weight <- matrix(NA,nyear, nspace,
                       dimnames = list(year = year, space = 1:nspace))
  Biomass.save <- matrix(NA,nyear, nspace,
                         dimnames= list(year = year, space = 1:nspace))
  Catch <- matrix(NA,nyear, dimnames = list(year = year))
  Catch.age <- matrix(NA,nage,nyear, dimnames = list(age = age, year = year))
  CatchN <- matrix(NA,nyear, dimnames = list(year = year))
  CatchN.age <- matrix(NA,nage,nyear, dimnames = list(age =age, year = year))


  R.save <- matrix(NA,nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  Fsel.save <- array(NA,dim = c(nage,nyear,nspace), dimnames = list(age = age, year = year, space = 1:nspace))
  Fseason.save <- array(NA,dim = c(nage, nyear, nspace,nseason), dimnames = list(age = age, year = year, space = 1:nspace,
                                                                                 season = 1:nseason))
  Fout.save <- array(NA, dim = c(nyear,nseason,nspace),
                     dimnames = list(year = year, season = 1:nseason, space = 1:nspace))

  N.save.age <- array(0,dim = c(nage,nyear+1, nspace, nseason),
                      dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  N.save.age.mid <- array(NA,dim = c(nage,nyear+1, nspace, nseason),
                          dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  R.save <- matrix(NA, nyear, nspace)
  V.save <- array(NA,dim = c(nyear, nspace, nseason), dimnames = list(
    year = year, space = 1:nspace, season = 1:nseason))

  E.save <- array(NA,dim = c(nyear, nspace, nseason), dimnames = list(
    year = year, space = 1:nspace, season = 1:nseason))


  Catch.save.age <- array(0,dim = c(nage,nyear, nspace, nseason),
                          dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))
  CatchN.save.age <- array(0,dim = c(nage,nyear, nspace, nseason),
                           dimnames = list(age = age, year = year, space = 1:nspace, season =1:nseason))

  survey <- array(NA,dim = c(nage, nyear, df$nsurvey))#,
                  #dimnames = list(year = year, age = df$age, survey = 1:df$nsurvey))

  survey.true <- array(NA, dim = c(nspace, nyear), dimnames = list(space = 1:nspace, year = year))
  surv.tot <- matrix(NA, nyear,nspace, dimnames = list(year = year, space = 1:nspace))

  age_comps_surv <- array(NA, dim = c(maxage,nyear), dimnames = list(age = 1:maxage,
                                                                            year = year)) #
  age_comps_surv_space <- array(NA, dim = c(maxage,nyear,nspace), dimnames = list(
    age = 1:maxage, year = year))

  N.survey <- matrix(NA,maxage ,nyear, dimnames = list(age = 1:maxage,
                                                              year= year))

  age_comps_catch <- array(NA, dim = c(maxage,nyear), dimnames = list(age = 1:maxage,
                                                                             year = year))
  age_comps_catch_space <- array(NA, dim = c(maxage,nyear,nspace), dimnames = list(
    age = 1:maxage, year = year, space = 1:nspace))

  age_comps_OM <- array(NA, dim = c(nage,nyear, nspace,nseason),
                        dimnames = list(age = age, year= year, space = 1:nspace, season = 1:nseason))

  Z.save <- array(NA, dim = c(df$nage, nyear,nspace,nseason), dimnames = list(age= age, year = year, space = 1:nspace,
                                                                              season = 1:nseason))

  Z.save[,1,1,1] <- M[,1,1]
  Catch.age[,1] <- 0 # Assumed no fishing before data started
  Catch[1]<- 0

  CatchN[1] <- 0
  CatchN.age[,1] <- 0

  survey[1] <- 1 # Surveys start later

  idx.save <- seq(1,tEnd, by = nseason)

  # Distribute over space
  if(is.null(df$Ninit)){
    Ninit <- N0
  }else{
    Ninit <- df$Ninit
  }
  #p.save <-matrix(NA,tEnd)


  for (space in 1:nspace){
    # if (season == 1){
    N.save.age[,1,space,1] <- Ninit # Just to initialize
    N.save.age.mid[,1,space,1] <- N.save.age[,1,space,1]*exp(-0.5*(M[,1,1]/nseason))
    # }else{
    #   N.save.age[,1,space,season] <- N.save.age[,1,space,season-1]*exp(-M/nseason)
    #   N.save.age.mid[,1,space,season] <- N.save.age[,1,space,season]*exp(-0.5*(M/nseason))
    # }
    # }
  }


  movemat <- array(0, dim = c(1,df$nage, df$nseason, nyear))


  for (yr in 1:(nyear)){ # Loop over years add one year for initial distribution

    w_catch <- df$weca[,yr,, drop = FALSE]
    w_catch[is.na(w_catch)] <- 0

    w_ssb <- df$west[,yr,, drop = FALSE]
    w_ssb[is.na(w_ssb)] <- 0


    sel <- df$sel[,yr,, drop = FALSE]
    sel[is.na(sel)] <- 0

    # Fyear <- F0[yr]*Fsel
    Myear <- M[,yr,, drop = FALSE] # Natural mortality
    mat.year <- mat[,yr,, drop = FALSE]


    for (season in 1:nseason){
      for (space in 1:nspace){

        if(season == 1){
         for(space in 1:nspace){
          SSB.weight[yr,space] <- sum(N.save.age[,yr,space,1]*w_ssb[,,season]*mat.year[,1,season], na.rm = TRUE)
          SSB[yr,space] <- SSB.weight[yr,space] #sum(N.save.age[,yr,space,1]*Mat.sel, na.rm = TRUE)

          SSB.all[1,space,1]<- sum(N.save.age[,1,space,1]*mat.year[,1,season], na.rm = TRUE)
        }

        }

        if(season == df$rseason){ # Recruitment season


            if(df$recruitment == 'Ricker'){

            R <- df$alpha*SSB[yr,]*exp(-betaSR*SSB[yr,])

            }



          if(df$recruitment == 'hockey'){
              R <- df$alpha + log(SSB[yr,space])

              if(SSB[yr, space] > df$betaSR){
                 R <- df$alpha+log(df$betaSR)
              }

              R <- exp(R)

              # add error
              Ry <- rnorm(1, mean = 0, sd = df$SDR)
              R.err <- exp(-0.5*df$b[yr]*SDR^2+Ry)

              N.save.age[1,yr,space,season] <- R*R.err
              R.save[yr,space] <- R*R.err

          }


        if(df$recruitment == 'estimated'){

          R <- df$Rin[yr]


          N.save.age[1,yr,space,season] <- R
          R.save[yr,space] <- R
          }

        }




        if(df$Fmodel == 'est'){
        Fseason <- F0[,yr, season]
        }

        if(df$Fmodel == 'sim'){

        Fseason <- df$Fin[yr, season]*sel[,1,season]

        }
        Mseason <- Myear[,1,season]


        Z <- Mseason+Fseason

        if(Z[1] == 0){
        Z[1] <- Z[2] # To avoid NaN calc
        }


        Z.save[,yr,space,season]<- Z
        Fseason.save[,yr,space, season] <- Fseason


        if(((space-1) == 0)){
          spaceidx <- 2
        }
        if(space == nspace){
          spaceidx <- nspace-1
        }
        if(space > 1 & space < nspace){
          spaceidx <- c(space-1,space+1)
        }

        if(df$move == FALSE){
          spaceidx <- 1
        }

        if(season <nseason){


          for(k in 1:length(spaceidx)){
            Nin.tmp<-  N.save.age[, yr,spaceidx[k],season]*exp(-Z)*(movemat[spaceidx[k],,season,yr])# add the ones come to the surrounding areas

            if(k == 1){
              Nin <- Nin.tmp
            }else{
              Nin <- Nin+Nin.tmp
            }

          }


          N.save.age[,yr,space,season+1] <- N.save.age[,yr,space,season]*exp(-Z)-
            N.save.age[, yr,space,season]*exp(-Z)*(movemat[space,,season,yr])+ # Remove the ones that leave
            Nin# add the ones come to the surrounding areas

          age_comps_OM[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])

          SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*mat.year[,1,season], na.rm = T)
          V.save[yr,space,season] <- sum(N.save.age[,yr,space,season]*sel[,1,season]*w_catch[,1,season])


          if(max(Fseason) > 0){
          Zcatch  <- Z
          Zcatch[1] <- Zcatch[2]

          Catch.save.age[, yr,space, season] <- (Fseason/(Zcatch))*(1-exp(-(Zcatch)))*N.save.age[,yr,space,season]*w_catch[,season]
          CatchN.save.age[, yr,space, season] <- (Fseason/(Zcatch))*(1-exp(-(Zcatch)))*N.save.age[,yr,space,season]


          E.save[yr,space,season] <- sum(Catch.save.age[,yr,space,season])/V.save[yr,space,season]


          }


        }else{

          for(k in 1:length(spaceidx)){
            Nin.tmp<-  + # Remove the ones that leave
              N.save.age[1:(nage-2), yr,spaceidx[k],season]*exp(-Z[1:(nage-2)])*(movemat[spaceidx[k],1:(nage-2),season,yr])## add the ones come to the surrounding areas

            Nin.plus.tmp <- (N.save.age[nage-1, yr,spaceidx[k],nseason]*exp(-Z[nage-1])+
                               N.save.age[nage, yr,spaceidx[k],nseason]*exp(-Z[nage]))*
              (movemat[spaceidx[k],nage, season,yr]) # Incoming



            if(k == 1){
              Nin <- Nin.tmp
              Nin.plus <- Nin.plus.tmp
            }else{
              Nin <- Nin+Nin.tmp
              Nin.plus <- Nin.plus.tmp+Nin.plus

            }

          }


          N.save.age[2:(nage-1),yr+1,space,1] <- N.save.age[1:(nage-2),yr,space,season]*exp(-Z[1:(nage-2)])-
            N.save.age[1:(nage-2), yr,space,season]*exp(-Z[1:(nage-2)])*(movemat[space,1:(nage-2),season,yr])+
            Nin #add the ones come to the surrounding areas
          # Plus group
          Nsurvive.plus <- (N.save.age[nage-1, yr,space, nseason]*exp(-Z[nage-1])+
                              N.save.age[nage, yr,space, nseason]*exp(-Z[nage]))

          Nout.plus <- Nsurvive.plus*(movemat[space,nage, season,yr]) # Leaving

          N.save.age[nage,yr+1,space,1] <- Nsurvive.plus- Nout.plus + Nin.plus


          age_comps_OM[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])

          SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*mat.year[,1,season], na.rm = T)

          V.save[yr,space,season] <- sum(N.save.age[,yr,space,season]*sel[,1,season]*w_catch[,1,season])


          if(max(Fseason) > 0){
          Zcatch  <- Z
          Zcatch[1] <- Zcatch[2]

          Catch.save.age[, yr,space, season] <- (Fseason/(Zcatch))*(1-exp(-(Zcatch)))*N.save.age[,yr,space,season]*w_catch[,1,season]
          CatchN.save.age[, yr,space, season] <- (Fseason/(Zcatch))*(1-exp(-(Zcatch)))*N.save.age[,yr,space,season]

          E.save[yr,space,season] <- sum(Catch.save.age[,yr,space,season])/V.save[yr,space,season]
          }

        }


      }


    } # End of season loop



    #Catch.age[,idx]  <- (Fyear/(Fyear+Myear))*(1-exp(-(Fyear+Myear)))*rowSums(N.save.age[,idx,,1])*w_catch # Calculate the catch in kg

    if(nseason>1){
      Catch.age[,yr] <- apply(Catch.save.age[,yr,,],MARGIN = 1,FUN = sum)
      Catch[yr] <- sum(Catch.save.age[,yr,,])

      CatchN.age[,yr] <- apply(CatchN.save.age[,yr,,],MARGIN = 1,FUN = sum)
      CatchN[yr] <- sum(CatchN.save.age[,yr,,])

    }else{

      if(nspace == 1){
        Catch.age[,yr] <- Catch.save.age[,yr,,]
        Catch[yr] <- sum(Catch.save.age[,yr,,])

        CatchN.age[,yr] <- CatchN.save.age[,yr,,]
        CatchN[yr] <- sum(CatchN.save.age[,yr,,])
      }else{
        Catch.age[,yr] <- rowSums(Catch.save.age[,yr,,])
        Catch[yr] <- sum(Catch.save.age[,yr,,])

        CatchN.age[,yr] <- rowSums(CatchN.save.age[,yr,,])
        CatchN[yr] <- sum(CatchN.save.age[,yr,,])
      }
    }

    for(surv in 1:df$nsurvey){

    if(df$move == FALSE){
      Nsurv <- N.save.age[,yr,,df$surveySeason]*
        exp(-.5*Z.save[,yr,space,df$surveySeason])*df$Q[,surv]
    }else{
      Nsurv <- rowSums(N.save.age[,yr,,df$surveySeason]*
                         exp(-Msurveymul*Z.save[,yr,space,df$surveySeason]))*df$Q[,surv]
    }


    survey[,yr,surv] <- N.save.age[,yr,,df$surveySeason[surv]]*
      exp(-.5*Z.save[,yr,space,df$surveySeason[surv]])*df$Q[,surv]*exp(rnorm(df$nage, mean =0, sd = 0.3)) # Change this


    survey[survey == 0] <- -1 # For TMB
    }



    Ntot.year <- Nsurv


    for(space in 1:nspace){
      Ntot.year <- N.save.age[,yr,space,df$surveyseason]

      if(nseason>1){
        Catch.tmp <- rowSums(CatchN.save.age[, yr,space,])
      }else{
        Catch.tmp <- CatchN.save.age[, yr,space,]
      }

      Catch.tot <- sum(CatchN.save.age[,yr,space,])

      age_comps_catch_space[1:(maxage-1),yr,space] <- Catch.tmp[2:(maxage)]/Catch.tot
      age_comps_catch_space[maxage,yr,space] <- sum(Catch.tmp[(maxage+1):nage])/Catch.tot



    }



  } # End of year loop
  #}

  if(df$move == FALSE){
    Nsave <- N.save.age[,,,nspace]
    SSB.save <- SSB
  }else{
    Nsave <- apply(N.save.age[,,,1],2,rowSums)
    SSB.save <- rowSums(SSB)
  }

  # Add names to output
  year_1 <- c(df$years,max(df$years+1))

  # Calulate Fbar


  #Fbar <- colSums(apply(Fseason.save[2:3,,1,], MARGIN = 2, colMeans))

  Fbar <- rep(0, nyear)

  for(time in 1:nyear){
        #if((age(i)>=df$Fbarage(0)) && (age(i)<=Fbarage(1))){

          Fbar[time] = sum(Fseason.save[(df$Fbarage[1]+1):(df$Fbarage[2]+1),time,,])/(df$Fbarage[2]-df$Fbarage[1]+1)

  }






  df.out   <- list(N.save = Nsave,
                   SSB = SSB,
                   N.save.age = N.save.age,
                   R.save = R.save,
                   V.save = V.save,
                   E.save = E.save,
                   SSB.all = SSB.all,
                   Catch.save.age = Catch.save.age,
                   CatchN.save.age = CatchN.save.age,
                   Catch = Catch,
                   Catch.age = Catch.age,
                   survey = survey,
                   Fbar = Fbar,
                   age_comps_OM = age_comps_OM,
                   age_catch = age_comps_catch,
                   Z = Z.save,
                   Fseason = Fseason.save,
                   Fsel = Fsel.save
               )

  return(df.out)

}
