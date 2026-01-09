#' Plot standard stock assessment figures
#'
#' @param OM smsR operating model
#' @param parms list of parameters
#' @param printFig TRUE FALSE return figure as list
#'
#' @return
#' Returns a figure with SSB, Fbar, recruitment and catch
#' @export
#'
#' @examples
#' smsPlots(df.tmb, sas)
#' smsPlots(df.tmb, sas, type = "stack", Blim = 110000)
#' smsPlots(df.tmb, sas, type = "wrap", Blim = 110000)
#' @importFrom ggplot2 ggplot aes geom_bar geom_line geom_hline geom_ribbon
#' @importFrom ggplot2 alpha scale_x_continuous scale_y_continuous coord_cartesian coord_flip
#' @importFrom ggplot2 facet_grid facet_wrap xlab ylab unit theme theme_classic theme_bw
#' @importFrom patchwork plot_layout
#'
plot.sms_om <- function(OM, printFig = 0) {


  if(ncol(OM$SSB) == 1){
  df.plot <- data.frame(
  SSB = as.numeric(OM$SSB),
  R = OM$R.save,
  Fbar = OM$Fbar,
  Catch = OM$Catch,
  years = as.numeric(rownames(OM$SSB))
  )
  }else{
    df.plot <- data.frame(
      SSB = apply(as.numeric(OM$SSB[,1]), 1, sum, na.rm =TRUE), # Sum over space
      R = OM$R.save,
      Fbar = OM$Fbar,
      Catch = OM$Catch,
      years = as.numeric(rownames(OM$SSB))
    )
  }
  # Plot SSB
  lims <- c(min(df.plot$SSB) / 2, max(df.plot$SSB) * 2) / 1000



    pssb <- ggplot(df.plot, aes(x = years, y = SSB / 1000)) +
      geom_line(linewidth = 1.2) +
      theme_classic() +
      scale_y_continuous("SSB\n(000 t)") +
      theme(legend.position.inside = c(0.8, .8)) +
      coord_cartesian(ylim = lims) +
      scale_x_continuous("")

    # Plot Recruitment
    # take care of crazy recruitment stuff
    lims <- c(min(df.plot$R) / 2, max(df.plot$R) * 2) / 1e6

    prec <- ggplot(df.plot, aes(x = years, y = R / 1e6)) +
      geom_line(linewidth = 1.2) +
      theme_classic() +
      scale_x_continuous("Year") +
      scale_y_continuous("Recruitment\n(millions)") +
      theme(legend.position = "none") +
      coord_cartesian(ylim = lims)


    # Plot Fishing mortality
    pF0 <- ggplot(df.plot, aes(x = years, y = Fbar)) +
      geom_line(size = 1.2) +
      theme_classic() +
      scale_x_continuous("Year") +
      scale_y_continuous("Fbar") +
      theme(legend.position = "none")


    # And catch
    pCatch <- ggplot(df.plot, aes(x = years, y = Catch / 1000)) +
      geom_line(size = 1.2) +
      theme_classic() +
      scale_x_continuous("") +
      scale_y_continuous("Catch\n(1000 t)") +
      theme(legend.position.inside = c(0.8, .8))
    # Plot SSB, recruitment, catch and fishing mortality
    pls <- (pssb + pCatch) / (prec + pF0)

    # pls
    if (printFig == 1) {
      print(pls)
    }


  #pls
  return(pls)
}
