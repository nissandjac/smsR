# getSSdata <- function(wd){
#
# #  wd <- "C:/Users/nsja/Documents/Github/test_smsR/petrale/petrale/models/2023.a034.001/"
#  require(r4ss)
#
#   dat <-  r4ss::SS_read(wd)
#
#   # Get catch data
#
#   Catch <- dat$dat$catch
#
#   # Get age comps
#   AC <- dat$dat$agecomp
#   ages <- dat$dat$agebin_vector # Get the min age fixed plz
#   nages <- length(ages)
#   years <- dat$dat$styr:dat$dat$endyr
#   nyears <- length(years)
#   nseason <- dat$dat$nseas
#
#   # Get weight at age in the catch
#   len <- dat$dat$lencomp
#   weight <- dat$dat$fleetinfo2
#
#   # turn catches into a catch at age
#   nfleets <- dat$dat$Nfleet
#   Catchobs_fleet <- array(NA, dim = c(nages, nyears, nseason, nfleets)) # No seasonality
#   weight_fleet <- array(NA, dim = c(nages, nyears, nseason, nfleets))
#
#   #
#
#
#
#
#   for(i in 1:unique(Catch$fleet)){
#     for(time in 1:nyears){
#       for(qrts in 1:nseason){
#
#       # Get the weight first
#       weight_
#
#
#       # Females only first. Add males later
#
#       if(years[time] < min(AC$Yr)){ # Take the last years comps
#       ACtmp <- AC %>% filter(Yr %in% min(AC$Yr))
#       }else{
#       ACtmp <- AC %>% filter(Yr %in% years[time])
#       }
#       # Get the index
#
#       nmstmp <- names(ACtmp)
#       idx <- nmstmp %in% paste('f',ages,sep='')
#       idx2 <- nmstmp %in% paste('m',ages,sep='')
#       # Get
#       Catchobs_fleet[,time,qrts,i] <- as.numeric(ACtmp[idx])*ACtmp$Nsamp
#
#       # Now multiply with the total catch
#
#       Catch.yr <- sum(Catch$catch[Catch$year%in% years[time]])
#
#
#       Catchobs_fleet[,time,qrts,i] <- Catchobs_fleet[,time,qrts,i] * Catch.yr
#
#       # Now divide by weight at age
#     }
#
#     }
#
#   }
#
#
#
