% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{getR}
\alias{getR}
\title{Estimated recruitment}
\usage{
getR(df.tmb, sas)
}
\arguments{
\item{df.tmb}{list of input parameters}

\item{sas}{fitted smsR model}
}
\value{
data frame of estimated recruitment.
}
\description{
Get a data frame of estimated recruitment from a fitted `smsR` model.
}
\details{
Last year is based on the Stock recruitment relationship chosen for the model fit. The estimated recruitment depends on the choice of recruitment season.
low and high are the 95\% confidence intervals.
SE is standard error of log R
}
\examples{
R <- getR(df.tmb, sas)
print(R)
}
