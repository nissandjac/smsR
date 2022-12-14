% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{getCatchN}
\alias{getCatchN}
\title{Get the number of catch individuals per age}
\usage{
getCatchN(df.tmb, sas)
}
\arguments{
\item{df.tmb}{input parameters}

\item{sas}{fitted smsR model}
}
\value{
data frame containing the numbers of individuals by age in the catch each year.
low and high are the 95\% confidence intervals.
SE is standard error of log N
}
\description{
Get the number of catch individuals per age
}
\examples{
CatchN <- getCatchN(df.tmb, sas)
print(CatchN)
}
