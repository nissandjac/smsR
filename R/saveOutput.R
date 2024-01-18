
#' Save output to evaluate assessment model
#'
#' @param df.tmb list of smsR input data
#' @param sas fitted smsR model
#' @param mr mohns rho object from the 'mohns_rho' function
#' @param name name of saved csv file
#' @param wd working directory
#' @param savefile TRUE or FALSE save the file to disc
#'
#' @return
#' Prints a csv file of diagnostic values
#' @export
#'
#' @examples
#'
#' dat <- saveOutput(df.tmb, sas, savefile = FALSE)
#' print(dat)
#'
#'
saveOutput <- function(df.tmb, sas, mr = NULL, savefile = TRUE,
                       name = 'summary', wd = getwd()){

  # Prepare a table with output stuff

  # Print s

  df.out <- getSummary(df.tmb, sas)

  if(savefile == TRUE){
    write.table(df.out, file = file.path(wd, paste(name,'.csv', sep = '')))
  }
  # add additional diagnostics
  ccv <- getCatchCV(df.tmb,sas)
  scv <- getSurveyCV(df.tmb, sas)

  # CV on hockey stick breakpoint??


  SDR<- as.numeric(sas$reps$par.fixed[names(sas$reps$par.fixed) == 'logSDrec'])


  model <- strsplit(wd, split = '/')[[1]]
  model <- model[length(model)]

  # Autocorrelation of residuals, age 1 and 2

  c <- getResidCatch(df.tmb, sas)
  c1 <- c[c$ages == 1 & c$season == 1,]

  lmmodel <- lm(ResidCatch~years, data = c1)
  ar_test <- lmtest::dwtest(lmmodel)

  c2 <- c[c$ages == 2 & c$season == 1,]
  lmmodel2 <- lm(ResidCatch~years, data = c2)
  ar_test2 <- lmtest::dwtest(lmmodel2)

  # Get the proper CV of SSB

  CV.ssb <- sas$reps$sd[names(sas$reps$value) == 'logSSB'][1:df.tmb$nyears]


  scv <- scv %>% dplyr::filter(surveyCV > 0)
  ccv <- ccv %>% dplyr::filter(catchCV > 0)


  df.indicators <- data.frame(mohns_r = as.numeric(MR$mohns[2]),
                   mohns_ssb = as.numeric(MR$mohns[1]),
                   mohns_F = as.numeric(MR$mohns[3]),
                   catchCV_1_avg = mean(ccv$catchCV[ccv$ages == 1]),
                   catchCV_2_avg = mean(ccv$catchCV[ccv$ages == 2]),
                   surveyCV_0 = mean(scv$surveyCV[scv$ages == 0]),
                   surveyCV_1 = mean(scv$surveyCV[scv$ages == 1]),
                   SDR = exp(SDR),
                   ARC1 = as.numeric(ar_test$statistic),
                   ARC2 = as.numeric(ar_test2$statistic),
                   model = model,
                   AIC = AIC(sas),
                   nparm = length(sas$reps$par.fixed),
                   SSB_sd_all = mean(CV.ssb),
                   SSB_sd_last3 = mean(CV.ssb[(df.tmb$nyears-2):df.tmb$nyears])

                   )
  row.names(df.indicators) <- NULL

  if(savefile == TRUE){
    write.table(df.indicators, file = file.path(wd, 'diagnostics.csv'), row.names = FALSE)
  }





return(df.indicators)
}
