#' Convert dataframes into matrices
#'
#' @param df data frame containing data
#' @param season number of quarters in model
#'
#' @return Returns the dataframe as a matrix for tmb
#' @export
#'
#' df <- canum_to_matrix(df, 1:4)
#'
#' @examples
#' df_to_matrix(df_lol, season = 1:4)
df_to_matrix <- function(df, season = 1:4){


  ages <- 0:max(df$Age)
  year <- unique(df$year)


  cmx <- array(NA, dim = c(length(ages), length(year), length(season)))

  #nn <- names(df)[4]


  for(i in 1:length(year)){
    for(j in 1:length(ages)){
      for(k in 1:length(season)){


        tmp <- df[df$Quarter == season[k] & df$year == year[i] & df$Age == ages[j],4]


        if(length(tmp)>0){
          cmx[j,i,k] <- tmp[[1]]

        }

      }
    }
  }



  return(cmx)
}
