
#' Organize TMB parameters in a list.
#'
#' @param mtrx matrix containing M, Mat, weca, west
#' @param Surveyobs matrix of survey observations
#' @param Catchobs matrix of catch observations
#' @param years vector of years
#' @param nseason number of seasons
#' @param nsurvey number of surveys
#' @param ages vector of ages
#' @param Fbarage vector of min and max ages used in Fbar calculations
#' @param recseason season where recruitment occurs
#' @param Fminage First age to start fishing
#' @param Fmaxage Last age with unique fishing mortality
#' @param Qminage vector length of nsurvey with minimum ages
#' @param Qmaxage vector length of nsurvey with max age with unique selectivity
#' @param Qlastage Last age observed in the survey
#' @param isFseason Season to calculate F in
#' @param CminageSeason  Minimum age with fishing mortality per season
#' @param endFseason  Last season last year with F
#' @param nocatch  Matrix size year x season with 1 or 0 showing catch
#' @param useEffort 1 or 0, use nominal CPUE to tune F
#' @param estimateCreep Estimate creep from nominal cpue
#' @param effort matrix size year x season that contains nominal CPUE
#' @param blocks Blocks with unique selectivity for fishing
#' @param surveyStart number between 0 and 1 determine the start of the survey within season
#' @param surveyEnd number between 0 and 1 determine the end of the survey within season
#' @param surveySeason vector determining which seasons surveys occur
#' @param minSDsurvey minimum allowed for all SDsurvey estimates
#' @param peneps epsilon used in penalty function for SDsurvey (try 1e-6 to 1e-3)
#' @param powers which ages are included in the power law calc
#' @param surveyCV list with ages with unique survey CVs. Length = nsurvey
#' @param catchCV list with ages with unique catch CVs. Length = nseason
#' @param recmodel chose recruitment model. 2 = hockeystick
#' @param estCV vector length 3 of which CVs to determine. 1-survey,2-catch, 3-stock recruitment
#' @param CVmin Currently not functional
#' @param betaSR breakpoint of SSB
#' @param nllfactor vector length 3 of weights of log likelihood. 1-survey, 2-catch,3-stock recruitment
#'
#' @return
#' list of stuff
#' @export
#'
#' @examples
get_TMB_parameters <- function(
  mtrx = NA,
  Surveyobs = NA,
  Catchobs = NA,
  years,
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
  endFseason = 4,
  nocatch = matrix(rep(1, nseason), nrow = length(years), ncol = nseason),
  useEffort = 0,
  estimateCreep = 0,
  effort = matrix(1, nrow = length(years), ncol = nseason),
  blocks = FALSE,
  surveyStart = rep(0, nsurvey),
  surveyEnd = rep(1, nsurvey),
  surveySeason = rep(1, nsurvey),
  minSDsurvey = 0.3,
  peneps = 1e-3,
  powers = list(NA),
  scv = array(0, dim = c(length(ages),length(years), nsurvey)),
  surveyCV = matrix(c(0,max(ages)), nrow = 2, ncol = nsurvey),
  catchCV = matrix(c(0,max(ages)), nrow = 2, ncol = nseason),
  recmodel = 2,
  estCV = c(0,0,0),
  CVmin = c(0.2,0.2,0.2),
  betaSR = NULL,
  nllfactor = rep(1,3)

){

  nsurvey <- dim(Surveyobs)[3]
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


  if(length(catchCV) > 1){
    for(i in 1:nseason){


      if(i == 1){
        no <- 1:length(catchCV[[i]])
      }else{
        no <- (max(no)+1):(max(no)+length(catchCV[[i]]))
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




      Cidx.CV[ages < CminageSeason[i],i] <- -98
      Cidx.CV[ages > CmaxageSeason[i],i] <- -98




    }

  }else{

    for(i in 1:nseason){


      no <- 1:length(catchCV[[i]])

      Cidx.CV[ages %in% catchCV[[i]],i] <- no
      Cidx.CV[ages > max(catchCV[[i]]),i] <- max(no)
      Cidx.CV[ages < min(catchCV[[i]]),i] <- min(no)

      Cidx.CV[ages < CminageSeason[i],nseason] <- -98
      Cidx.CV[ages > CmaxageSeason[i],nseason] <- -98

    }




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

#  Catchobs[Catchobs <= 1] <- 0





  df.tmb <- list(
    weca = mtrx$weca,
    west = mtrx$west,
    Surveyobs = Surveyobs,
    Catchobs = Catchobs,
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
    effort = effort.in,
    bidx = bidx,
    useBlocks = useBlocks,
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
    M = mtrx$M,
    Mat = mtrx$mat,
    scv = scv,
    surveyStart = surveyStart,
    surveyEnd = surveyEnd,#c(0.1,1,0.001),
    surveySeason = surveySeason,
    minSDsurvey = minSDsurvey,
    peneps = peneps,
    powers = powersexp,
    recmodel = recmodel, # 1 is hockey stick
    estCV = estCV,
    CVmin = CVmin,
    betaSR = betaSR,
    nllfactor = nllfactor

  )

  return(df.tmb)

}
