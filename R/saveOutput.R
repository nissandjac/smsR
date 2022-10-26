
#' Save output to evaluate assessment model
#'
#' @param df.tmb
#' @param sas
#' @param MR
#' @param save
#' @param Fbarage
#' @param name
#' @param wd
#'
#' @return
#' @export
#'
#' @examples
saveOutput <- function(df.tmb, sas, MR = NULL, save = TRUE,
                      Fbarage = c(1,2), name = 'summary', wd = getwd()){

  # Prepare a table with output stuff

  # Print s

  df.out <- getSummary(df.tmb, sas, Fbarage)

  if(save == TRUE){
    write.table(df.out, file = file.path(wd, paste(name,'.csv', sep = '')))
  }
  # add additional diagnostics
  ccv <- getCatchCV(df.tmb,sas)
  scv <- getSurveyCV(sas)

  # CV on hockey stick breakpoint??


  SDR<- sas$reps$par.fixed[names(sas$reps$par.fixed) == 'logSDrec']


  model <- strsplit(wd, split = '/')[[1]]
  model <- model[length(model)]

  # Autocorrelation of residuals, age 1 and 2

  c <- getResidCatch(df.tmb, sas)
  c1 <- c1[c1$ages == 1,]

  ac1 <- acf(c1$ResidCatch)


  c2 <- c[c$ages == 2,]

  ac2 <- acf(c1$ResidCatch)




  df.indicators <- data.frame(mohns_r = as.numeric(MR$mohns[2]),
                   mohns_ssb = as.numeric(MR$mohns[1]),
                   mohns_F = as.numeric(MR$mohns[3]),
                   catchCV_1_avg = mean(ccv$catchCV[ccv$ages == 1]),
                   catchCV_2_avg = mean(ccv$catchCV[ccv$ages == 2]),
                   surveyCV_0 = mean(scv$surveyCV[scv$ages == 0]),
                   surveyCV_1 = mean(scv$surveyCV[scv$ages == 1]),
                   SDR = exp(SDR),
                   model = model
                   )


  if(save == TRUE){
    write.table(df.indicators, file = file.path(wd, 'diagnostics.csv'), row.names = FALSE)
  }





return(df.indicators)
}
