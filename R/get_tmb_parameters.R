
#' Prepare data input for an smsR model
#'
#' @param mtrx matrix containing weca, west, mat, and M
#' @param Surveyobs Survey observations
#' @param Catchobs Catch observations
#' @param propM proportion M before spawning
#' @param propF proportion F before spawning
#' @param years years with available data
#' @param startYear start year of assessment (optional)
#' @param endYear end year of assessment (optional)
#' @param nseason number of seasons
#' @param nsurvey number of surveys
#' @param ages age span
#' @param Fbarage Ages to average F to Fbar
#' @param recseason recruitment season
#' @param Fminage minimum age to calculate F
#' @param Fmaxage Maximum age to calculate F
#' @param Qminage Minimum age included in survey
#' @param Qmaxage Max age included in survey
#' @param Qlastage Last age with unique survey selectivity
#' @param isFseason seasons to calculate F
#' @param CminageSeason minimum age to use catchdata
#' @param CmaxageSeason Maximum age to use catchdata
#' @param endFseason Last season with fishing
#' @param nocatch input matrix of 1 or 0 defining if F should be calculated
#' @param useEffort Use nominal effort to calculate F
#' @param estimateCreep Estimate technological creep as a parameter. Needs useEffort = 1
#' @param effort matrix of nominal effort per season
#' @param blocks unique selectivity blocks
#' @param surveyStart fraction into the season where the survey starts
#' @param surveyEnd fraction into the season where the survey ends
#' @param surveySeason season where survey occurs
#' @param leavesurveyout vector of 1 and 0 to exclude surveys
#' @param minSDsurvey minimum CV of survey
#' @param minSDcatch minimum CV of catch
#' @param peneps parameter for regulating minimum CV for survey
#' @param penepsC parameter for regulating minimum CV of catch
#' @param powers apply powerlaw for density dependent survey observations
#' @param scv add time varying survey CV (input matrix)
#' @param surveyCV survey cv grouping
#' @param catchCV catch cv grouping
#' @param recmodel recruitment model
#' @param estCV which CVs to estimate as parameters
#' @param CVmin minimum CVs
#' @param betaSR hockey stick break point
#' @param nllfactor negative log likelihood weighting
#' @param randomF try random effect fishing mortality
#'
#'
#' @return
#' returns a
#' @export
#'
#' @examples
#' get_TMB_paramters(mtrx = sandeel_1r$lhs,
#'                  Surveyobs = sandeel_1r$survey,
#'                  Catchobs = sandeel_1r$Catch ,
#'                  years = 1983:2022) # Not run
#'
get_TMB_parameters <- function(

  mtrx = NULL,
  Surveyobs = NULL,
  Catchobs = NULL,
  propM = NULL,
  propF = NULL,
  years,
  startYear = min(years),
  endYear = max(years),
  nseason = 4,
  nsurvey = 2,
  ages = 0:20,
  Fbarage = c(1,max(ages)),
  recseason = 1,
  Fminage = 0,
  Fmaxage = max(ages),
  Qminage = rep(0, length(ages)),
  Qmaxage = rep(max(ages), length(ages)),
  Qlastage = Qmaxage,
  isFseason = rep(1, nseason),
  CminageSeason = rep(0, nseason),
  CmaxageSeason = rep(max(ages), nseason),
  endFseason = nseason,
  nocatch = matrix(rep(1, nseason), nrow = length(years), ncol = nseason),
  useEffort = 0,
  estimateCreep = 0,
  effort = matrix(1, nrow = length(years), ncol = nseason),
  blocks = FALSE,
  surveyStart = rep(0, nsurvey),
  surveyEnd = rep(1, nsurvey),
  surveySeason = rep(1, nsurvey),
  leavesurveyout = rep(1,nsurvey),
  minSDsurvey = 0.3,
  minSDcatch = 0.2,
  peneps = 1e-3,
  penepsC = 1e-3,
  powers = list(NA),
  scv = array(0, dim = c(length(ages),length(years), nsurvey)),
  surveyCV = matrix(c(0,max(ages)), nrow = 2, ncol = nsurvey),
  catchCV = matrix(c(0,max(ages)), nrow = 2, ncol = nseason),
  recmodel = 2,
  estCV = c(0,0,0),
  CVmin = c(0.2,0.2,0.2),
  betaSR = NULL,
  nllfactor = rep(1,3),
  randomF = 0

){


  # Remove surveys for sensitivity analysis
    if(sum(leavesurveyout) != nsurvey){




      Surveyobs <- Surveyobs[,,leavesurveyout == 1, drop = FALSE]
      Qminage <- Qminage[leavesurveyout == 1]
      Qmaxage <- Qmaxage[leavesurveyout == 1]
      Qlastage <- Qlastage[leavesurveyout == 1]
      powers <- powers[leavesurveyout == 1]
      surveyCV <- surveyCV[leavesurveyout == 1]

    }





  nsurvey <- dim(Surveyobs)[3]
  if(nsurvey == 0){
    warning('probably doesnt work without survey')
  }



  Qidx <- rep(0, nsurvey)

  if(nsurvey > 1){
  for(i in 2:nsurvey){
    Qidx[i] <- length(ages[ages == Qminage[i-1]]:ages[ages == Qlastage[i-1]])+Qidx[i-1]

  }}else{
    Qidx <- 0
  }



  # Fix the survey CV groups
  nage <- length(ages)
  Qidx.CV <- matrix(0, nage, nsurvey)
  no <- 1
  maxage <- max(ages)

  for(k in 1:nsurvey){

      tmpCV <- surveyCV[[k]]+1 # Go from age to index
      vec <- rep(0, nage)

      if(length(tmpCV) == 1){

        vec[tmpCV:(maxage+1)] <- no
        no <- no+1

      }else{



      for(i in 1:(length(tmpCV))){

        if(i < length(tmpCV)){
          tmp.idx <- tmpCV[i]:(tmpCV[i+1]-1)
          vec[tmp.idx] <- no
          no <- no+1
        }else{
          vec[tmpCV[length(tmpCV)]:(tmpCV[i]+1)] <- no
          no <- no+1
        }
      }

        vec[max(tmpCV):nage] <- vec[max(tmpCV)]

      }

      # Expand the tmp CV into a nage length  vector

      rm.idx <- which(0:maxage < Qminage[k] | 0:maxage > Qmaxage[k])
      vec[rm.idx] <- -98
      Qidx.CV[,k] <- vec
  }


  Qidx.CV <- Qidx.CV-1 # Scale to c++ indexing







  nyear <- length(years)
  # Turn the block into an index
  effort.in <- matrix(0, nyear, nseason)


  if(blocks[1] != FALSE){
    blocks.idx <- which(years %in% blocks) #
    # Extend the index to years
    bidx <- rep(NA, nyear)

    for(i in 1:(length(blocks.idx)-1)){

      len <- blocks.idx[i]:(blocks.idx[i+1]-1)
      bidx[len] <- i-1


    }
    bidx[blocks.idx[length(blocks.idx)]:nyear] <- length(blocks.idx)-1


    # Change mean effort to 1 (within blocks)
    nblocks <- length(blocks.idx)

    for(i in 1:nblocks){
      tmpidx <- which((bidx+1) == i)

      tmpeffort <- effort[tmpidx,]

      Meffort <- sum(tmpeffort)/length(tmpeffort[tmpeffort>0])

      effort.in[tmpidx,] <- effort[tmpidx, ]/Meffort

      # }
    }

  }else{

    tmpeffort <- effort
    Meffort <- mean(tmpeffort[tmpeffort > 0])
    effort.in<- effort/Meffort
    bidx <- rep(0, nyear)




  }

  # Scale effort to 1

  # Add index to determine if b has several blocks

  if(all(blocks == FALSE)){
    useBlocks <- 0
  }else{
    useBlocks <- 1
  }

  # if(all(effort != 1)){
  #   useEffort <- 1
  # }

  if(nseason > 1){
  isFseason[length(isFseason)] <- 0 # This is a weird standard thing in sms
  }
  # Do the power law calcs
  if(is.na(powers[[1]])){
    powersexp <- matrix(0, nrow = length(ages), ncol = nsurvey)
  }else{
    powersexp <- matrix(0, nrow = length(ages), ncol = nsurvey)

    for(i in 1:nsurvey){
      if(is.na(powers[[i]]) == 0){
      powersexp[powers[[i]]+1,i] <-1
      }
    }
  }



  # Fix the survey CV groups
  Cidx.CV <- matrix(NA, nage, nseason)

  if(min(CminageSeason) > Fminage){
    Fminage <- min(CminageSeason)
  }


  #if(length(catchCV) > 1){
    for(i in 1:nseason){

      if(min(catchCV[[i]]) != 0){
        catchCV[[i]] <- c(0,catchCV[[i]])
      }


      if(i == 1){

#        if(min(catchCV[[i]])>0){
        no <- 1:length(catchCV[[i]])
        # }else{
        # no <- 0:(length(catchCV[[i]])-1)
        # }
        }else{
          no <- (max(no)+1):(max(no)+length(catchCV[[i]]))

          no <- no - CminageSeason[i]


        }
      Cidx.CV[ages %in% catchCV[[i]],i] <- no
      Cidx.CV[ages > max(catchCV[[i]]),i] <- max(no)
      Cidx.CV[ages < min(catchCV[[i]]),i] <- min(no)


      # Do a loop for NA check
      for(a in 2:nage){
        if(is.na(Cidx.CV[a, i])){
          Cidx.CV[a,i] <- Cidx.CV[a-1,i]
        }
      }

      if(CminageSeason[i] < Fminage){
        CminageSeason[i] <- Fminage
      }


      Cidx.CV[ages < CminageSeason[i],i] <- -98
      Cidx.CV[ages > CmaxageSeason[i],i] <- -98




    }
#
#   }else{
#
#     for(i in 1:nseason){
#
#
#       no <- 1:length(catchCV[[i]])
#
#       Cidx.CV[ages %in% catchCV[[i]],i] <- no
#       Cidx.CV[ages > max(catchCV[[i]]),i] <- max(no)
#       Cidx.CV[ages < min(catchCV[[i]]),i] <- min(no)
#
#       for(a in 2:nage){
#         if(is.na(Cidx.CV[a, i])){
#           Cidx.CV[a,i] <- Cidx.CV[a-1,i]
#         }
#       }
#
#       Cidx.CV[ages < CminageSeason[i],nseason] <- -98
#       Cidx.CV[ages > CmaxageSeason[i],nseason] <- -98
#
#     }
#
#
#
#
#   }

  if(min(Cidx.CV[Cidx.CV>0]) == 1){ # Do some index fixing
    Cidx.CV <- Cidx.CV + 1
  }

  Cidx.CV <- Cidx.CV - 2 # Convert to C++ idx

  CVgroups <- NA

  for(i in 1:length(catchCV)){
    CVgroups[i] <- length(catchCV[[i]])

  }


  # Do Catch CV for internal CV calcs
  # CCV.out <- array(unlist(catchCV), dim = c())

  if(length(unique(CVgroups))> 1 & estCV[2] == 2){
    stop('sms currently doesnt support multiple catch CVs between seasons with internally calculated SD')
  }

  catchCVout <- matrix(rep(catchCV[[1]], nseason), ncol= nseason)
#
#   if(estCV[2] == 2){
#   catchCVout <- catchCV
#
#   }

  # Calc the number of catch observations

  no <- matrix(0, nrow = nrow(catchCVout), ncol = nseason)

  for(i in 1:nrow(catchCVout)){
    for(qrts in 1:nseason){
      if((i-1) >= CminageSeason[qrts]){

        if(i < nrow(catchCVout)){
        idx <- catchCVout[i]:(catchCVout[i+1]-1)+1
        }else{
        idx <- (catchCVout[i]+1):nage
        }

        Out <- Catchobs[idx,,qrts]

        no[i,qrts] <- length(Out[Out >0])
      }
    }
  }

  if(nrow(catchCVout) == 1 & nseason == 1){
    no[i,qrts] <- length(Catchobs[Catchobs > 0])
  }

#  Catchobs[Catchobs <= 1] <- 0


  if(is.null(propM)){
    propM <- array(0, dim = c( nage , nyear+1, nseason))
  }

  if(is.null(propF)){
    propF <- array(0, dim = c( nage , nyear+1, nseason))
  }


  if(startYear > min(years)){
    weca <- mtrx$weca[,c(years %in% startYear:max(years),TRUE),]
    west <- mtrx$west[,c(years %in% startYear:max(years),TRUE),]
    M <- mtrx$M[,c(years %in% startYear:max(years),TRUE),]
    Mat <- mtrx$mat[,c(years %in% startYear:max(years),TRUE),]
    propM <- propM[,c(years %in% startYear:max(years),TRUE),]
    propF <- propF[,c(years %in% startYear:max(years),TRUE),]

    Surveyobs <- Surveyobs[,which(years %in% startYear:max(years)),]
    Catchobs <- Catchobs[,which(years %in% startYear:max(years)),]

    scv <- scv[,which(years %in% startYear:max(years)),]
    effort <- effort.in[which(years %in% startYear:max(years)),]
    nocatch <- nocatch[which(years %in% startYear:max(years)),]
    bidx <- bidx[which(years %in% startYear:max(years))]

    years <- startYear:max(years)
    nyears <- length(years)

  }else{
    weca <- mtrx$weca
    west <- mtrx$west
    M <- mtrx$M
    Mat <- mtrx$mat
    effort <- effort.in
  }


  if(endYear < max(years)){
    weca <- mtrx$weca[,c(years %in% startYear:endYear,TRUE),]
    west <- mtrx$west[,c(years %in% startYear:endYear,TRUE),]
    M <- mtrx$M[,c(years %in% startYear:endYear,TRUE),]
    Mat <- mtrx$mat[,c(years %in% startYear:endYear,TRUE),]
    propM <- propM[,c(years %in% startYear:endYear,TRUE),]
    propF <- propF[,c(years %in% startYear:endYear,TRUE),]

    Surveyobs <- Surveyobs[,which(years %in% startYear:endYear),]
    Catchobs <- Catchobs[,which(years %in% startYear:endYear),]

    scv <- scv[,which(years %in% startYear:endYear),]
    effort <- effort.in[which(years %in% startYear:endYear),]
    nocatch <- nocatch[which(years %in% startYear:endYear),]
    bidx <- bidx[which(years %in% startYear:endYear)]

    years <- startYear:endYear
    nyears <- length(years)

  }



  df.tmb <- list(
    weca = weca,
    west = west,
    Surveyobs = Surveyobs,
    Catchobs = Catchobs,
    propM = propM,
    propF = propF,
    no = no,
    years = years,
    age = ages,
    Fbarage = Fbarage,
    nage = length(ages),
    nseason = nseason,
    nyears = length(years),
    nsurvey = dim(Surveyobs)[3],
    recseason = recseason,
    useEffort = useEffort,
    estimateCreep = estimateCreep,
    effort = effort,
    bidx = bidx,
    useBlocks = useBlocks,
    blocks = blocks,
    Fminage = Fminage,
    Fmaxage = Fmaxage,
    Qminage = Qminage,
    Qmaxage = Qmaxage,
    Qlastage = Qlastage,
    Qidx = Qidx,
    Qidx_CV = Qidx.CV,
    Cidx_CV = Cidx.CV,
    catchCV = catchCVout,
    isFseason = isFseason, # Fishing mortality in how many quarterS? ,
    endFseason = endFseason,
    CminageSeason = CminageSeason,
    nocatch = nocatch,
    M = M,
    Mat = Mat,
    scv = scv,
    surveyStart = surveyStart,
    surveyEnd = surveyEnd,#c(0.1,1,0.001),
    surveySeason = surveySeason,
    minSDsurvey = minSDsurvey,
    minSDcatch = minSDcatch,
    surveyCV = surveyCV,
    catchCVin = catchCV,
    peneps = peneps,
    penepsC = penepsC,
    powers = powersexp,
    recmodel = recmodel, # 1 is hockey stick
    estCV = estCV,
    CVmin = CVmin,
    betaSR = betaSR,
    nllfactor = nllfactor,
    randomF = randomF

  )

  return(df.tmb)

}
