#' Title
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
getDataSandeel <- function(wd,
                           maxage = NA,
                           survey.age = list(),
                           survey.years = list(),
                           survey.names = NA,
                           survey.quarter = NA,
                           effort = 0,
                           years,
                           seasons){


  nyears <- length(years)

  canum <- read.table(file.path(wd,'canum.in'),
                      comment.char = "#", skip = 0, header = FALSE)

  colnames(canum) <- 0:maxage
  canum$year <- rep(years, each = length(seasons))
  canum$Quarter <- rep(seasons, nyears)
  canum <- canum %>% pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'catchage')


  survey <- read.table(file.path(wd, 'fleet_catch.in'), fill = TRUE, col.names = c('effort',paste(1:4, sep = ',')))

  ages <- 0:maxage



  survey.out <- matrix(NA, nrow = nrow(survey), ncol = length(ages) + 4) # Plus efficiency, name, quarter and years
  survey.out <- as.data.frame(survey.out)
  colnames(survey.out) <- c('eff', ages, 'Survey','Quarter','year')
  survey.out$eff <- survey[,1]


  for(k in 1:length(survey.names)){
    sname <-  rep(survey.names[k], length(survey.years[[k]]))
    syears <- survey.years[[k]]


    if(k == 1){
      idx <- 1:length(syears)
    }else{
      idx <- (idx[length(idx)]+1):((idx[length(idx)])+length(syears))
    }


    survey.out$Survey[idx] <- sname
    survey.out$year[idx] <- syears
    survey.out$Quarter[idx] <- rep(survey.quarter[k], length(idx))


    for(i in 1:length(survey.age[[k]])){
      survey.out[idx, names(survey.out) == survey.age[[k]][i]] <- survey[idx, i + 1]

    }



  }

  survey.out <- survey.out %>% tidyr::pivot_longer(2:(length(ages)+1), names_to = 'Age', values_to = 'cpue')

  survey <- survey.out
  survey$cpue[is.na(survey$cpue)] <- -1

  if(effort == TRUE){
    effort <- read.table(file.path(wd, 'effort.in'), header = FALSE, fill =TRUE)
    colnames(effort) <- 'effort'
    effort$year <- rep(years, length(seasons))
    effort$Quarter <- rep(seasons, each = length(years))
  }else{effort =NA}

  # SCale the effort to 1
  # # Change mean effort to 1 (within blocks)
  # effort.in <- matrix(0, df.tmb$nyears, df.tmb$nseason)
  #
  # for(i in 1:nblocks){
  #   tmpidx <- which(df.tmb$bidx == i)
  #
  #   tmpeffort <- effort[tmpidx,]
  #
  #   Meffort <- sum(tmpeffort)/length(tmpeffort[tmpeffort>0])
  #
  #   effort.in[tmpidx,] <- effort[tmpidx, ]/Meffort
  #
  #   # }
  # }
  #


  # Get the life history parameters

  #
  M <- read.table(file.path(wd,'natmor.in'),
                            comment.char = "#", skip = 0, header = FALSE)

  colnames(M) <- 0:maxage

  # Remove the last couple of years (it goes to 2023)
  M <- M[1:(nyears*max(seasons)),]

  M$year <- rep(years, each = length(seasons))
  M$Quarter <- rep(seasons, nyears)
  M <- M %>% pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'mortality')

  # Maturity

  mat <- read.table(file.path(wd,'propmat.in'),
                         comment.char = "#", skip = 0, header = FALSE)

  colnames(mat) <- 0:maxage

  # Remove the last couple of years (it goes to 2023)
  mat <- mat[1:(nyears*max(seasons)),]

  mat$year <- rep(years, each = length(seasons))
  mat$Quarter <- rep(seasons, nyears)
  mat <- mat %>% pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'mortality')

  # Now west and weca
  weca <- read.table(file.path(wd,'weca.in'),
                  comment.char = "#", skip = 0, header = FALSE)

  colnames(weca) <- 0:maxage

  # Remove the last couple of years (it goes to 2023)
  weca <- weca[1:(nyears*max(seasons)),]

  weca$year <- rep(years, each = length(seasons))
  weca$Quarter <- rep(seasons, nyears)
  weca <- weca %>% pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'weca')


  west <- read.table(file.path(wd,'west.in'),
                     comment.char = "#", skip = 0, header = FALSE)

  colnames(west) <- 0:maxage

  # Remove the last couple of years (it goes to 2023)
  west <- west[1:(nyears*max(seasons)),]

  west$year <- rep(years, each = length(seasons))
  west$Quarter <- rep(seasons, nyears)
  west <- west %>% pivot_longer(1:(maxage+1), names_to = 'Age', values_to = 'west')


  # Make them into matrices

  mtrx <- list(M = df_to_matrix(M, season = seasons),
               mat = df_to_matrix(mat, season = seasons),
               west = df_to_matrix(west, season = seasons),
               weca = df_to_matrix(weca, season= seasons))





  return(list(canum = canum,
              survey = survey,
              effort = effort,
              mtrx = mtrx)
  )
}
