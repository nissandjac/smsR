#' Turn data frames into multidimensional matrices for TMB
#'
#' @param df.in data frame with info
#' @param season number of seasons to standardize over
#'
#' @return
#' returns a matrix
#' @export
#'
#' @examples
#'

df_to_matrix_multi <- function(df.in, season = 1:4){
  
  
  ages <- 0:max(df.in$Age)
  year <- unique(df.in$year)
  spp <- unique(df.in$spp)
  
  
  cmx <- array(NA, dim = c(length(ages), 
                           length(year), 
                           length(season), 
                           length(spp)))
  
  
  
  for(i in 1:length(year)){
    for(j in 1:length(ages)){
      for(k in 1:length(season)){
        for(l in 1:length(spp)){
        tmp <- df.in[df.in$Quarter == season[k] & df.in$year == year[i] & df.in$Age == ages[j] & df.in$spp == spp[l], 5]
        
        if(length(tmp)>0){
          cmx[j,i,k,l] <- tmp[[1]]
          
        }
        
        }
      }
    }
  }
  
  
  return(cmx)
}
