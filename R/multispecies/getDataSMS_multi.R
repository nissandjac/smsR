#' Read data from an SMS model
#'
#' @param wd working directory
#' @param maxage maxage
#' @param survey.age list showing ages in surveys
#' @param survey.years years ssurveys occur
#' @param survey.names names of survey
#' @param survey.quarter which quarters do surveys happen
#' @param effort set to true to read in effort data
#' @param years years to calculate
#' @param seasons seasons to calc
#'
#' @return
#' @export
#'
#' @examples
#' x <- getDataSandeel # Not run
getDataSMS_multi <- function(wd, spp = NULL, is.dynamic = 1,
                       maxage = NA,
                       survey.age = list(),
                       survey.years = list(),
                       survey.names = NA,
                       survey.quarter = NA,
                       effort = FALSE,
                       years,
                       seasons,
                       scv.tv = 0){
  
  
  nyears <- length(years)
  nseason <- length(seasons)
  spp_m <- sum(is.dynamic)
  nage <- length(0:maxage)
  
  canum <- read.table(file.path(wd,'canum.in'),
                      comment.char = "#", skip = 0, header = FALSE, nrows = nyears * nseason * spp_m)
  #canum <- canum[-nrow(canum), ]
  
  colnames(canum) <- 0:maxage
  canum$year <- rep(rep(years, each = length(seasons)), spp_m)
  canum$Quarter <- rep(seasons, nyears*spp_m)
  canum$spp <- rep(spp[which(is.dynamic == 1)], each = nyears * nseason)
  
  
  canum <- canum %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'catchage') %>% mutate(Age = as.numeric(as.character(Age)))
  
  
  # Try to read  the fleet info.dat 
  
  finfo <- read.table(file.path(wd, 'fleet_info.dat'), fill = TRUE, flush = TRUE, col.names = 1:10)
  
  nsurveys <- finfo[2:(spp_m+1),1]
  
  eidx <- spp_m+1
  
  survey.years <- survey.start <- survey.age <- list()
  power.model <- survey.season <- survey.nvar <- survey.lastage <- rep(NA, length(spp))
  survey.names <- janitor::make_clean_names(rep(spp[which(is.dynamic == 1)], nsurveys))
  # survey.quarter = NA,
  # 
  
  for(i in 1:sum(nsurveys)){
    
    survey.years[[i]] <- finfo[eidx+i,1]:finfo[eidx+i,2]
    survey.start[[i]] <- c(finfo[eidx+i,3],finfo[eidx+i,4])
    survey.age[[i]] <- c(finfo[eidx+i,5],finfo[eidx+i,6])
    survey.lastage[i] <- finfo[eidx+i,7]
    power.model[i] <- finfo[eidx+i,8]
    survey.season[i] <-finfo[eidx+i,9]
    survey.nvar[i] <-finfo[eidx+i,10]
  }
  
  
  survey <- read.table(file.path(wd, 'fleet_catch.in'), fill = TRUE, col.names = c('effort',paste(1:max(unlist(survey.age)), sep = ',')))
  syears <- c(years, max(years)+1) # Some 2020 stuff going on here 
  
  for(i in 1:spp_m){
    
    for(j in 1:nsurveys[i]){
    
      if(i == 1 & j == 1){
      val <- 1
    }else{
      val <- val+1
    } 
      
    if(i == 1 & j == 1){
      idx <- 1
    }
      
      
    tmp <- data.frame(years = rep(syears, each = nage),  
                      age = rep(0:maxage, nyears+1),
                      CPUE = -1,
                      effort = -1,
                      survey = survey.names[val]
                      )
    
    
    yr.tmp <- survey.years[[val]]
    age.tmp <- survey.age[[val]][1]:survey.age[[val]][2]
    
    cpue.tmp <- as.matrix(survey[idx:(idx+length(yr.tmp)-1),2:(length(age.tmp)+1)], ncol = length(age.tmp))
    eff.tmp <- survey[idx:(idx+length(yr.tmp)-1),1]
    
    
    for(k in 1:ncol(cpue.tmp)){
      tmp$CPUE[tmp$age == age.tmp[k] & tmp$years %in% yr.tmp] <- cpue.tmp[,k]
      tmp$effort[tmp$age == age.tmp[k] & tmp$years %in% yr.tmp] <- eff.tmp
      
    }
    
    idx <- (idx+length(yr.tmp))
    
    
    if(i == 1 & j == 1){
      survey.out <- tmp
    }else{
      survey.out <- rbind(survey.out, tmp)
    }
    
    }
  }
  
  survey.list <- list(CPUE = survey.out ,
                      survey.season = survey.season,
                      survey.names = survey.names, 
                      power.model = power.model, 
                      survey.nvar = survey.nvar,
                      survey.start = survey.start,
                      lastage = survey.lastage)
  # Cut the bottom if there's only one survey
 
  
  
  # Get the life history parameters
  
  #
  M <- read.table(file.path(wd,'natmor.in'),
                  comment.char = "#", skip = 0, header = FALSE, nrows = nyears * nseason * spp_m)
  
  colnames(M) <- 0:maxage
  
  # Remove the last couple of years
  M$year <- rep(rep(years, each = length(seasons)), spp_m)
  M$Quarter <- rep(seasons, nyears*spp_m)
  M$spp <- rep(spp[which(is.dynamic == 1)], each = nyears * nseason)
  
  M <- M %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'mortality') %>% mutate(Age = as.numeric(as.character(Age)))
  
  # Maturity
  
  mat <- read.table(file.path(wd,'propmat.in'),
                    comment.char = "#", skip = 0, header = FALSE, nrows = nyears * nseason * spp_m)
  
  colnames(mat) <- 0:maxage
  
  # Remove the last couple of years (it goes to 2023)
  
  # Remove the last couple of years
  mat$year <- rep(rep(years, each = length(seasons)), spp_m)
  mat$Quarter <- rep(seasons, nyears*spp_m)
  mat$spp <- rep(spp[which(is.dynamic == 1)], each = nyears * nseason)
  
  mat <- mat %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'maturity') %>% mutate(Age = as.numeric(as.character(Age)))
  
  
  weca <- read.table(file.path(wd,'weca.in'),
                    comment.char = "#", skip = 0, header = FALSE, nrows = nyears * nseason * spp_m)
  
  colnames(weca) <- 0:maxage
  
  # Remove the last couple of years (it goes to 2023)
  
  # Remove the last couple of years
  weca$year <- rep(rep(years, each = length(seasons)), spp_m)
  weca$Quarter <- rep(seasons, nyears*spp_m)
  weca$spp <- rep(spp[which(is.dynamic == 1)], each = nyears * nseason)
  
  weca <- weca %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'weca') %>% mutate(Age = as.numeric(as.character(Age)))
  
  
  
  west <- read.table(file.path(wd,'west.in'),
                     comment.char = "#", skip = 0, header = FALSE, nrows = nyears * length(spp)*nseason)
  
  colnames(west) <- 0:maxage
  
  # Remove the last couple of years (it goes to 2023)
  
  # Remove the last couple of years
  west$year <- rep(rep(years, each = length(seasons)), length(spp))
  west$Quarter <- rep(seasons, nyears*length(spp))
  west$spp <- rep(spp, each = nyears * nseason)
  
  west <- west %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'west') %>% mutate(Age = as.numeric(as.character(Age)))
  
  
  
  # Make them into matrices
  
  mtrx <- list(M = df_to_matrix_multi(M, season = seasons),
               mat = df_to_matrix_multi(mat, season = seasons),
               west = df_to_matrix_multi(west, season = seasons),
               weca = df_to_matrix_multi(weca, season= seasons))
  
  
  # Load the stomach content info 
  
  # Other predators N at age 
  
  other_pred <- spp[which(is.dynamic == 0)]
  
  op_N <- read.table(file.path(wd,'other_pred_N.in'),
                     comment.char = "#", skip = 0, header = FALSE, nrows = nyears * length(other_pred)*nseason)
  
  
  colnames(op_N) <- 0:maxage
  op_N$year <- rep(rep(years, each = length(seasons)), length(other_pred))
  op_N$Quarter <- rep(seasons, nyears*length(other_pred))
  op_N$spp <- rep(other_pred, each = nyears * nseason)
  
  op_N <- op_N %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'op_N') %>% mutate(Age = as.numeric(as.character(Age)))
  
  # Get some relevant species info (might export later)
  spp.info <- read.table(file.path(wd, 'SMS.dat'), , fill = TRUE, flush = TRUE, col.names = 1:length(spp))
  is.predator <- spp.info[15:(15+length(spp)-1),7]
  
  other_food <- read.table(file.path(wd, 'other_food.in'), nrow = 1)
  names(other_food) <- spp[which(is.predator > 0)]
  
  
  
  # Mean length of predators 
  west$year <- rep(rep(years, each = length(seasons)), length(spp))
  west$Quarter <- rep(seasons, nyears*length(spp))
  west$spp <- rep(spp, each = nyears * nseason)
  
  west <- west %>% tidyr::pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'west') %>% mutate(Age = as.numeric(as.character(Age)))
  
  
  
  
  return(list(canum = canum,
              survey = survey.list,
              mtrx = mtrx)
  )
}
