#' Title
#'
#' @param df.tmb list of tmb parameters
#' @param sas estimated model.
#' @param type default, stack, or wrap
#' @param Blim optional value to graph if type is stack or wrap
#'
#' @return
#' @export
#'
#' @examples
#' smsPlots(df.tmb, sas)
#' smsPlots(df.tmb, sas, type="stack", Blim=110000)
#' smsPlots(df.tmb, sas, type="wrap", Blim=110000)
#' @importFrom ggplot2 ggplot aes geom_bar geom_line geom_hline geom_ribbon
#' @importFrom ggplot2 alpha scale_x_continuous scale_y_continuous coord_cartesian coord_flip
#' @importFrom ggplot2 facet_grid facet_wrap xlab ylab unit theme theme_classic theme_bw
#'
 smsPlots <- function(df.tmb, sas, type="default", Blim=NULL, printFig = TRUE){

  SSB <- getSSB(df.tmb, sas)
  rec <- getR(df.tmb,sas)
  Fbar <- getFbar(df.tmb, sas)
  Catch <- getCatch(df.tmb,sas)

  if(type=="default"){
    # Plot SSB
    pssb <- ggplot(SSB, aes(x = years, y = SSB))+geom_line(size = 1.4)+
        theme_classic()+geom_ribbon(aes(ymin = low, ymax = high), fill = alpha('red', 0.2), linetype = 0)+
        scale_y_continuous('SSB')+theme(legend.position = c(0.8,.8))

    # Plot Recruitment
    # take care of crazy recruitment stuff
    lims <- c(min(rec$R)/2, max(rec$R)*2)

    prec <- ggplot(rec, aes(x = years, y = R))+geom_line(size = 1.4)+
        theme_classic()+geom_ribbon(aes(ymin = low, ymax = high), fill = alpha('red', 0.2), linetype = 0)+
        scale_y_continuous('recruitment')+theme(legend.position = 'none')+coord_cartesian(ylim = lims)


    # Plot Fishing mortality
    pF0 <-  ggplot(Fbar, aes(x = years, y = Fbar))+geom_line(size = 1.3)+theme_classic()+
        geom_ribbon(aes(ymin = low, ymax = high), fill = alpha('red', 0.2), linetype = 0)+
        scale_y_continuous('fishing mortality')+theme(legend.position = 'none')


    # And catch
    pCatch <-  ggplot(Catch, aes(x = years, y = Catch))+geom_line(size = 1.4)+
        theme_classic()+geom_ribbon(aes(ymin = low, ymax = high), fill = alpha('red', 0.2), linetype = 0)+
        scale_y_continuous('Catch')+theme(legend.position = c(0.8,.8))
    # Plot SSB, recruitment, catch and fishing mortality
    pls <- (pssb + pCatch)/(prec + pF0)

    #pls
    if(printFig == 1){
      print(pls)
    }

    return(pls)
  }#end default plot

  if(type %in% c("stack", "wrap")){

    #Make names to attach them
    SSB$variable <- "SSB \n(1000 tonnes)"
    rec$variable <- "Recruitment \n(billion individuals)"
    Fbar$variable <- paste0("F[ages ", df.tmb$Fbarage[1], "-", df.tmb$Fbarage[2], "]")
    Catch$variable <- "Catch \n(1000 tonnes)"

    names(SSB)[1] <- "mid"
    names(rec)[1] <- "mid"
    names(Fbar)[1] <- "mid"
    names(Catch)[1] <- "mid"

    #Scale to match names and stock assessment sheet
    SSB$mid <- SSB$mid/1000
    rec$mid <- rec$mid/1000000
    Catch$mid <- Catch$mid/1000
    SSB$low <- SSB$low/1000
    rec$low <- rec$low/1000000
    Catch$low <- Catch$low/1000
    SSB$high <- SSB$high/1000
    rec$high <- rec$high/1000000
    Catch$high <- Catch$high/1000


    results <- rbind(Fbar, SSB, rec, Catch)

    #FIXME: don't repeat yourself, find smarter way to deal with missing Blim
    if(!is.null(Blim)){
      Blimdat <- data.frame(years=df.tmb$years,
                          Blim=Blim/1000,
                          lo=NA, hi=NA,
                          variable="SSB \n(1000 tonnes)", type=NA)
      if(type=="stack"){
        pls <- ggplot(results, aes(x=years, y=mid))+
          geom_line()+geom_ribbon(aes(ymin=low, ymax=high), alpha=.4, linetype=0)+
          facet_grid(variable~., scale="free", switch="both")+ylab(NULL)+
          theme(panel.spacing = unit(0, "lines"))+
          geom_hline(data=Blimdat, aes(yintercept=Blim), lty=2)+
          #	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
          theme(legend.position="none")+
          scale_x_continuous(expand = c(0,0))+xlab(NULL)+
          theme_bw()
      }
      if(type=="wrap"){
        pls <- ggplot(results, aes(x=years, y=mid))+
          geom_line()+geom_ribbon(aes(ymin=low, ymax=high), alpha=.4, linetype=0)+
          facet_wrap(~variable, scale="free")+ylab(NULL)+
          theme(panel.spacing = unit(0, "lines"))+
          geom_hline(data=Blimdat, aes(yintercept=Blim), lty=2)+
          #	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
          theme(legend.position="none")+
          scale_x_continuous(expand = c(0,0))+xlab(NULL)+
          theme_bw()
      }
    }#end Blim
    if(is.null(Blim)){
      if(type=="stack"){
        pls <- ggplot(results, aes(x=years, y=mid))+
          geom_line()+geom_ribbon(aes(ymin=low, ymax=high), alpha=.4, linetype=0)+
          facet_grid(variable~., scale="free", switch="both")+ylab(NULL)+
          theme(panel.spacing = unit(0, "lines"))+
          #	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
          theme(legend.position="none")+
          scale_x_continuous(expand = c(0,0))+xlab(NULL)+
          theme_bw()
      }
      if(type=="wrap"){
        pls <- ggplot(results, aes(x=years, y=mid))+
          geom_line()+geom_ribbon(aes(ymin=low, ymax=high), alpha=.4, linetype=0)+
          facet_wrap(~variable, scale="free")+ylab(NULL)+
          theme(panel.spacing = unit(0, "lines"))+
          #	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
          theme(legend.position="none")+
          scale_x_continuous(expand = c(0,0))+xlab(NULL)+
          theme_bw()
      }
    }
  }
  pls
  return(pls)
}
