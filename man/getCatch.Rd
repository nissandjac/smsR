% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{getCatch}
\alias{getCatch}
\title{Estimated catch}
\usage{
getCatch(df.tmb, sas)
}
\arguments{
\item{df.tmb}{input parameters from \code{\link{get_TMB_parameters}}}

\item{sas}{fitted smsR object from \code{\link{runAssessment}}}
}
\value{
data frame containing the biomass of catch each year.
}
\description{
Retrieve the total estimated catch weight from a fitted `smsR` object. Standard unit is kg.
}
\details{
low and high are the 95\% confidence intervals.
SE is standard error of log Catch
}
\examples{
catch <- getCatch(df.tmb, sas)
print(catch)
}
\seealso{
\code{\link{getYield}}, \code{\link{getCatchN}}
}
