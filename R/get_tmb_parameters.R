#' Prepare Input Data for an smsR Model
#' @description
#'  This function prepares the necessary input data for an `smsR` model.
#'  Required input are `mtrx`, `Surveyobs`,`Catchobs`, and `year`.
#'  Note that the dimensions of the input data needs to be consistent with year, season and ages.
#'
#' @param mtrx Matrix containing life history parameters: `weca`, `west`, `mat`, and `M`.
#' @param Surveyobs Matrix of survey observations (numbers at age).
#' @param Catchobs Matrix of catch observations (numbers at age).
#' @param propM Proportion of natural mortality (`M`) occurring before spawning.
#' @param propF Proportion of fishing mortality (`F`) occurring before spawning.
#' @param years Vector of years with available data.
#' @param startYear Optional. Start year of the assessment.
#' @param endYear Optional. End year of the assessment.
#' @param nseason Number of seasons in the model.
#' @param nsurvey Number of surveys included in the model.
#' @param ages Vector specifying the age span of the data.
#' @param Fbarage Vector of ages used to calculate the mean fishing mortality (`Fbar`).
#' @param recseason Season in which recruitment occurs.
#' @param Fminage Minimum age used to calculate fishing mortality (`F`).
#' @param Fmaxage Maximum age used to calculate fishing mortality (`F`).
#' @param Qminage Minimum age included in the survey data.
#' @param Qmaxage Maximum age included in the survey data.
#' @param Qlastage Last age with unique survey selectivity.
#' @param isFseason Vector specifying seasons in which fishing mortality (`F`) is calculated.
#' @param CminageSeason Minimum age for using catch data.
#' @param CmaxageSeason Maximum age for using catch data.
#' @param endFseason Last season in which fishing occurs.
#' @param nocatch Matrix (1/0) defining whether fishing mortality (`F`) should be calculated.
#' @param useEffort Logical. If `TRUE`, nominal effort is used to calculate `F`.
#' @param estimateCreep Logical. If `TRUE`, estimates technological creep as a parameter (requires `useEffort = TRUE`).
#' @param effort Matrix of nominal fishing effort per season.
#' @param blocks Vector of unique selectivity blocks.
#' @param surveyStart Fraction into the season when the survey starts.
#' @param surveyEnd Fraction into the season when the survey ends.
#' @param surveySeason Season in which the survey occurs.
#' @param leavesurveyout Vector (1/0) indicating which surveys to exclude.
#' @param minSDsurvey Minimum coefficient of variation (CV) for survey data.
#' @param minSDcatch Minimum CV for catch data.
#' @param maxSDcatch Maximum CV for catch data.
#' @param peneps Parameter regulating the minimum CV for surveys.
#' @param penepsC Parameter regulating the minimum CV for catch data.
#' @param penepsCmax Parameter regulating the maximum CV for catch data.
#' @param powers Logical. If `TRUE`, applies a power-law function for density-dependent survey observations.
#' @param scv Matrix specifying time-varying survey CV.
#' @param surveySD Grouping variable for survey SD.
#' @param catchSD Grouping variable for catch SD.
#' @param recmodel Recruitment model specification. Options:
#'   - `1`: Hockey-stick model with input breakpoint (`betaSR`).
#'   - `2`: Recruitment estimated as deviations from an estimated mean value (can include environmental input).
#'   - `3`: Beverton-Holt model.
#' @param estSD Vector indicating which deviation parameters to estimate.
#' @param CVmin Minimum CV values for estimation.
#' @param betaSR Hockey-stick model breakpoint parameter.
#' @param nllfactor Weighting factor for the negative log-likelihood.
#' @param randomF Logical. If `TRUE`, fishing mortality (`F`) is treated as a random effect.
#' @param randomM Logical. If `TRUE`, natural mortality (`M`) is estimated as a time-varying random effect.
#' @param randomR Logical. If `TRUE`, recruitment deviations are estimated as a random effect.
#' @param nenv Number of environmental covariates affecting recruitment.
#' @param Mprior Prior standard deviation (SD) for `M` deviations from the first year of the time series.
#' @param M_min Minimum age included in the time-varying `M` estimation.
#' @param M_max Maximum age included in the time-varying `M` estimation.
#' @param MCV Age distribution of the CV for time-varying `M`.
#' @param SDMprior Prior SD for `M` CV estimation.
#' @param Pred_in Matrix of predator inputs for MICE modeling.
#'
#' @details
#'  The required inputs include:
#'   - `mtrx`: A matrix containing life history parameters (i.e,, `weca`, `west`, `mat`, `M`).
#'   - `Catchobs`: A matrix of observed catch numbers at age.
#'   - `Surveyobs`: A matrix of observed survey numbers at age.
#'
#' Additionally, the user must specify the `years` corresponding to the input data.

#'
#'
#' @return A list containing all required input data for the smsR model.
#' @export
#'
#' @examples
#' get_TMB_parameters(
#'   mtrx = sandeel_1r$lhs,
#'   Surveyobs = sandeel_1r$survey,
#'   Catchobs = sandeel_1r$Catch,
#'   years = 1983:2022
#' ) # Not run

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
    Fbarage = c(1, max(ages)),
    recseason = 1,
    Fminage = 0,
    Fmaxage = max(ages),
    Qminage = rep(0, nsurvey),
    Qmaxage = rep(max(ages), nsurvey),
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
    leavesurveyout = rep(1, nsurvey),
    minSDsurvey = 0.3,
    minSDcatch = 0.2,
    maxSDcatch = sqrt(2),
    peneps = 1e-3,
    penepsC = 1e-3,
    penepsCmax = 1e-3,
    Mprior = 0.5,
    SDMprior = 0.2,
    powers = list(NA),
    Pred_in = NA,
    M_min = 0,
    M_max = max(ages),
    scv = array(0, dim = c(length(ages), length(years), nsurvey)),
    surveySD = matrix(c(0, max(ages)), nrow = 2, ncol = nsurvey),
    catchSD = matrix(c(0, max(ages)), nrow = 2, ncol = nseason),
    MCV = matrix(c(0, max(ages)), nrow = 2, ncol = 1),
    estSD = c(0, 0, 0),
    CVmin = c(0.2, 0.2, 0.2),
    betaSR = 0,
    nllfactor = rep(1, 3),
    randomF = 0,
    randomM = 0,
    randomR = 0,
    recmodel = 1,
    nenv = 0,
    warn = FALSE) {
  # Remove surveys for sensitivity analysis

   if (sum(leavesurveyout) != nsurvey) {
    Surveyobs <- Surveyobs[, , leavesurveyout == 1, drop = FALSE]
    Qminage <- Qminage[leavesurveyout == 1]
    Qmaxage <- Qmaxage[leavesurveyout == 1]
  #  Qlastage <- Qlastage[leavesurveyout == 1]
    powers <- powers[leavesurveyout == 1]
    surveySD <- surveySD[leavesurveyout == 1]
    nsurvey <- sum(leavesurveyout)
    surveySeason <- surveySeason[leavesurveyout == 1]
  }


  if(length(Pred_in) == 1){
    nalphaM <- length(MCV[[1]])
  }else{
    nalphaM <- nrow(Pred_in)
  }

  if(randomM == 1 & all(is.na(Pred_in))){
    nalphaM_CV <- 1
  }else{
    nalphaM_CV <- nrow(Pred_in)
  }

  nsurvey <- dim(Surveyobs)[3]
  if (nsurvey == 0) {
    warning("probably doesnt work without survey")
  }


  for(i in 1:nsurvey){
    if(surveyEnd[i] == surveyStart[i]){
      surveyEnd[i] <- surveyStart[i]*1.001 # Small hack

    }
  }


  Qidx <- rep(0, nsurvey)

  if (nsurvey > 1) {
    for (i in 2:nsurvey) {
      Qidx[i] <- length(ages[ages == Qminage[i - 1]]:ages[ages == Qlastage[i - 1]]) + Qidx[i - 1]
    }
  } else {
    Qidx <- 0
  }



  # Fix the survey CV groups
  nage <- length(ages)
  Qidx.CV <- matrix(0, nage, nsurvey)
  no <- 1
  maxage <- max(ages)

  for (k in 1:nsurvey) {
    tmpCV <- surveySD[[k]] + 1 # Go from age to index
    vec <- rep(0, nage)

    if (length(tmpCV) == 1) {
      vec[tmpCV:(maxage + 1)] <- no
      no <- no + 1
    } else {
      for (i in 1:(length(tmpCV))) {
        if (i < length(tmpCV)) {
          tmp.idx <- tmpCV[i]:(tmpCV[i + 1] - 1)
          vec[tmp.idx] <- no
          no <- no + 1
        } else {
          vec[tmpCV[length(tmpCV)]:(tmpCV[i] + 1)] <- no
          no <- no + 1
        }
      }

      vec[max(tmpCV):nage] <- vec[max(tmpCV)]
    }

    # Expand the tmp CV into a nage length  vector

    rm.idx <- which(0:maxage < Qminage[k] | 0:maxage > Qmaxage[k])
    vec[rm.idx] <- -98
    Qidx.CV[, k] <- vec
  }


  Qidx.CV <- Qidx.CV - 1 # Scale to c++ indexing


  # Fix the survey CV groups

  if(randomM == 1){

    Midx.CV <- matrix(0, nage, nalphaM_CV)
    no <- 1
    maxage <- max(ages)


    for (k in 1:nalphaM_CV) {
    tmpCV <- MCV[[k]] + 1 # Go from age to index
    vec <- rep(0, nage)

    if (length(tmpCV) == 1) {
      vec[tmpCV:(maxage + 1)] <- no
      no <- no + 1
    } else {
      for (i in 1:(length(tmpCV))) {
        if (i < length(tmpCV)) {
          tmp.idx <- tmpCV[i]:(tmpCV[i + 1] - 1)
          vec[tmp.idx] <- no
          no <- no + 1
        } else {
          vec[tmpCV[length(tmpCV)]:(tmpCV[i] + 1)] <- no
          no <- no + 1
        }
      }

      vec[max(tmpCV):nage] <- vec[max(tmpCV)]
    }

    # Expand the tmp CV into a nage length  vector

    rm.idx <- which(0:maxage < M_min[k] | 0:maxage > M_max[k])
    vec[rm.idx] <- -98
    Midx.CV[, k] <- vec
  }


  Midx.CV <- Midx.CV - 1 # Scale to c++ indexing

  M_nparms <- length(unique(Midx.CV[Midx.CV > -1]))

  }else{
    Midx.CV <- matrix(1, 1)
    M_nparms <- 1
  }

  nyear <- length(years)
  # Turn the block into an index
  effort.in <- matrix(0, nyear, nseason)


  if (blocks[1] != FALSE) {
    blocks.idx <- which(years %in% blocks) #
    # Extend the index to years
    bidx <- rep(NA, nyear)

    for (i in 1:(length(blocks.idx) - 1)) {
      len <- blocks.idx[i]:(blocks.idx[i + 1] - 1)
      bidx[len] <- i - 1
    }
    bidx[blocks.idx[length(blocks.idx)]:nyear] <- length(blocks.idx) - 1


    # Change mean effort to 1 (within blocks)
    nblocks <- length(blocks.idx)

    for (i in 1:nblocks) {
      tmpidx <- which((bidx + 1) == i)

      tmpeffort <- effort[tmpidx, ]

      Meffort <- sum(tmpeffort) / length(tmpeffort[tmpeffort > 0])

      effort.in[tmpidx, ] <- effort[tmpidx, ] / Meffort

      # }
    }
  } else {
    tmpeffort <- effort
    Meffort <- mean(tmpeffort[tmpeffort > 0])
    effort.in <- effort / Meffort
    bidx <- rep(0, nyear)
  }

  # Scale effort to 1

  # Add index to determine if b has several blocks

  if (all(blocks == FALSE)) {
    useBlocks <- 0
  } else {
    useBlocks <- 1
  }

  # if(all(effort != 1)){
  #   useEffort <- 1
  # }

  if (nseason > 1) {
    isFseason[length(isFseason)] <- 0 # This is a weird standard thing in sms
  }
  # Do the power law calcs
  if (is.na(powers[[1]])) {
    powersexp <- matrix(0, nrow = length(ages), ncol = nsurvey)
  } else {
    powersexp <- matrix(0, nrow = length(ages), ncol = nsurvey)

    for (i in 1:nsurvey) {
      if (is.na(powers[[i]]) == 0) {
        powersexp[powers[[i]] + 1, i] <- 1
      }
    }
  }



  # Fix the survey CV groups
  Cidx.CV <- matrix(NA, nage, nseason)

  if (min(CminageSeason) > Fminage) {
    Fminage <- min(CminageSeason)
  }


  # if(length(catchSD) > 1){
  for (i in 1:nseason) {
    if (min(catchSD[[i]]) != 0) {
      catchSD[[i]] <- c(0, catchSD[[i]])
    }


    if (i == 1) {
      #        if(min(catchSD[[i]])>0){
      no <- 1:length(catchSD[[i]])
      # }else{
      # no <- 0:(length(catchSD[[i]])-1)
      # }
    } else {
      no <- (max(no) + 1):(max(no) + length(catchSD[[i]]))

      no <- no - CminageSeason[i]
    }
    Cidx.CV[ages %in% catchSD[[i]], i] <- no
    Cidx.CV[ages > max(catchSD[[i]]), i] <- max(no)
    Cidx.CV[ages < min(catchSD[[i]]), i] <- min(no)


    # Do a loop for NA check
    for (a in 2:nage) {
      if (is.na(Cidx.CV[a, i])) {
        Cidx.CV[a, i] <- Cidx.CV[a - 1, i]
      }
    }

    if (CminageSeason[i] < Fminage) {
      CminageSeason[i] <- Fminage
    }


    Cidx.CV[ages < CminageSeason[i], i] <- -98
    Cidx.CV[ages > CmaxageSeason[i], i] <- -98
  }
  #
  #   }else{
  #
  #     for(i in 1:nseason){
  #
  #
  #       no <- 1:length(catchSD[[i]])
  #
  #       Cidx.CV[ages %in% catchSD[[i]],i] <- no
  #       Cidx.CV[ages > max(catchSD[[i]]),i] <- max(no)
  #       Cidx.CV[ages < min(catchSD[[i]]),i] <- min(no)
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

  if (min(Cidx.CV[Cidx.CV > 0]) == 1) { # Do some index fixing
    Cidx.CV <- Cidx.CV + 1
  }

  Cidx.CV <- Cidx.CV - 2 # Convert to C++ idx

  CVgroups <- NA

  for (i in 1:length(catchSD)) {
    CVgroups[i] <- length(catchSD[[i]])
  }


  # Do Catch CV for internal CV calcs
  # CCV.out <- array(unlist(catchSD), dim = c())

  if (length(unique(CVgroups)) > 1 & estSD[2] == 2) {
    stop("sms currently doesnt support multiple catch CVs between seasons with internally calculated SD")
  }

  catchSDout <- matrix(rep(catchSD[[1]], nseason), ncol = nseason)
  #
  #   if(estSD[2] == 2){
  #   catchSDout <- catchSD
  #
  #   }

  # Calc the number of catch observations

  no <- matrix(0, nrow = nrow(catchSDout), ncol = nseason)

  for (i in 1:nrow(catchSDout)) {
    for (qrts in 1:nseason) {
      if ((i - 1) >= CminageSeason[qrts]) {
        if (i < nrow(catchSDout)) {
          idx <- catchSDout[i]:(catchSDout[i + 1] - 1) + 1
        } else {
          idx <- (catchSDout[i] + 1):nage
        }

        Out <- Catchobs[idx, , qrts]

        no[i, qrts] <- length(Out[Out > 0])
      }
    }
  }

  if (nrow(catchSDout) == 1 & nseason == 1) {
    no[i, qrts] <- length(Catchobs[Catchobs > 0])
  }

  #  Catchobs[Catchobs <= 1] <- 0


  if (is.null(propM)) {
    propM <- array(0, dim = c(nage, nyear + 1, nseason))
  }

  if (is.null(propF)) {
    propF <- array(0, dim = c(nage, nyear + 1, nseason))
  }

  if(is.null(mtrx$M)){
    stop('Provide natural mortality')
  }


  if(is.null(dim(mtrx$M))){
    mtrx$M <- matrix(rep(mtrx$M, each = nage), nrow = nage, ncol = nyear+1)
    message('Turning M vector into age matrix')
  }


  if(length(mtrx$mat) == nage){
    mtrx$mat <- matrix(mtrx$mat, nrow = nage, ncol = nyear)
  }

  # Add a projection year to weca and west

  if(dim(mtrx$weca)[2] == nyear){

    weca.mean <- apply(mtrx$weca[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    mtrx$weca <- abind::abind(mtrx$weca, array(weca.mean, dim = c(nage, 1, nseason) ), along = 2)

    if(warn == TRUE){
    warning('add projection year to weca. Using average of last 3 years')
    }

  }


  if(dim(mtrx$west)[2] == nyear){

    west.mean <- apply(mtrx$west[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    mtrx$west <- abind::abind(mtrx$west,array(west.mean, dim = c(nage, 1, nseason) ), along = 2)
    if(warn == TRUE){

    warning('add projection year to weca. Using average of last 3 years')

    }
  }


  if(dim(mtrx$mat)[2] == nyear){

    mat.mean <- apply(mtrx$mat[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    mtrx$mat <- abind::abind(mtrx$mat,array(mat.mean, dim = c(nage, 1, nseason) ), along = 2)
    if(warn == TRUE){

    warning('add projection year to mat. Using average of last 3 years')
}
  }

  if(dim(mtrx$M)[2] == nyear){

    M.mean <- apply(mtrx$M[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    mtrx$M <- abind::abind(mtrx$M,array(M.mean, dim = c(nage, 1, nseason) ), along = 2)
    if(warn == TRUE){

    warning('add projection year to M. Using average of last 3 years')

    }
  }

  if(dim(propM)[2] == nyear){

    propM.mean <- apply(propM[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    propM <- abind::abind(propM,array(M.mean, dim = c(nage, 1, nseason) ), along = 2)
    if(warn == TRUE){

    warning('add projection year to propM. Using average of last 3 years')

    }
  }

  if(dim(propF)[2] == nyear){

    propF.mean <- apply(propF[,(nyear-2):nyear,, drop =FALSE], FUN = rowMeans, MAR = c(3))
    propF <- abind::abind(propF,array(M.mean, dim = c(nage, 1, nseason) ), along = 2)
    if(warn == TRUE){

    warning('add projection year to propF. Using average of last 3 years')

    }
  }

  # Add dimension names to variables

  dnames <- list(ages = ages, years = c(years, max(years)+1), seasons = 1:nseason)

  dimnames(mtrx$weca) <- dnames
  dimnames(mtrx$west) <- dnames
  dimnames(mtrx$M) <- dnames
  dimnames(mtrx$mat) <- dnames
  dimnames(propM) <- dnames
  dimnames(propF) <- dnames

  dnames <- list(ages = ages, years = c(years), seasons = 1:nseason)
  dimnames(Catchobs) <- dnames

  if(is.null(dimnames(Surveyobs))){
  dnames <- list(ages = ages, years = c(years), seasons = 1:nsurvey)
  dimnames(Surveyobs) <- dnames
  }
  # Some other things
  dimnames(effort.in) <- list(years = years, seasons = 1:nseason)
  dimnames(nocatch) <- list(years = years, seasons = 1:nseason)


  if (startYear > min(years)) {
    weca <- mtrx$weca[, c(years %in% startYear:max(years), TRUE), ]
    west <- mtrx$west[, c(years %in% startYear:max(years), TRUE), ]
    M <- mtrx$M[, c(years %in% startYear:max(years), TRUE), ]
    Mat <- mtrx$mat[, c(years %in% startYear:max(years), TRUE), ]
    propM <- propM[, c(years %in% startYear:max(years), TRUE), ]
    propF <- propF[, c(years %in% startYear:max(years), TRUE), ]

    Surveyobs <- Surveyobs[, which(years %in% startYear:max(years)), ]
    Catchobs <- Catchobs[, which(years %in% startYear:max(years)), ]

    scv <- scv[, which(years %in% startYear:max(years)), ]
    effort <- effort.in[which(years %in% startYear:max(years)), ]
    nocatch <- nocatch[which(years %in% startYear:max(years)), ]
    bidx <- bidx[which(years %in% startYear:max(years))]

    years <- startYear:max(years)
    nyears <- length(years)
  } else {
    weca <- mtrx$weca
    west <- mtrx$west
    M <- mtrx$M
    Mat <- mtrx$mat
    effort <- effort.in
  }


  if (endYear < max(years)) {
    weca <- mtrx$weca[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]
    west <- mtrx$west[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]
    M <- mtrx$M[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]
    Mat <- mtrx$mat[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]
    propM <- propM[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]
    propF <- propF[, c(years %in% startYear:endYear, TRUE), , drop = FALSE]

    Surveyobs <- Surveyobs[, which(years %in% startYear:endYear), , drop = FALSE]
    Catchobs <- Catchobs[, which(years %in% startYear:endYear), , drop = FALSE]

    scv <- scv[, which(years %in% startYear:endYear), , drop = FALSE]
    effort <- effort.in[which(years %in% startYear:endYear), , drop = FALSE]
    nocatch <- nocatch[which(years %in% startYear:endYear), , drop = FALSE]
    bidx <- bidx[which(years %in% startYear:endYear), drop = FALSE]

    years <- startYear:endYear
    nyears <- length(years)
  }


 # Fill environmental matrix with 0's if its not there
  env_matrix <- matrix(0,  nenv,length(years))
  M_matrix <- matrix(0, nalphaM, length(years))



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
    catchSD = catchSDout,
    Midx_CV = Midx.CV,
    isFseason = isFseason, # Fishing mortality in how many quarterS? ,
    endFseason = endFseason,
    CminageSeason = CminageSeason,
    nocatch = nocatch,
    M = M,
    Pred_in = Pred_in,
    Mat = Mat,
    scv = scv,
    surveyStart = surveyStart,
    surveyEnd = surveyEnd, # c(0.1,1,0.001),
    surveySeason = surveySeason,
    minSDsurvey = minSDsurvey,
    minSDcatch = minSDcatch,
    maxSDcatch = maxSDcatch,
    catchSDin = catchSD,
    surveySD = surveySD,
    peneps = peneps,
    penepsC = penepsC,
    penepsCmax = penepsCmax,
    Mprior = Mprior,
    SDMprior = SDMprior,
    nalphaM = nalphaM,
    M_nparms = M_nparms,
    powers = powersexp,
    recmodel = recmodel, # 1 is hockey stick
    estSD = estSD,
    CVmin = CVmin,
    betaSR = betaSR,
    nllfactor = nllfactor,
    randomF = randomF,
    randomR = randomR,
    randomM = randomM,
    nenv = nenv,
    env_matrix = env_matrix,# Number of environmental parameters,
    M_matrix = M_matrix
  )

  # Check if everything looks right
  check_tmb_error(df.tmb)

  return(df.tmb)
}
