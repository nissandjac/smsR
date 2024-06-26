% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{getM}
\alias{getM}
\title{get the natural mortality at age from the input data and plot it}
\usage{
getM(df.tmb, plotFig = FALSE)
}
\arguments{
\item{df.tmb}{smsR list of input parameters}

\item{plotFig}{plot the data (TRUE OR FALSE)}
}
\value{
returns a data frame of weight at aeg
}
\description{
get the natural mortality at age from the input data and plot it
}
\examples{

M <- getM(df.tmb, plotFig = TRUE) # Plot the natural mortality by age

}
