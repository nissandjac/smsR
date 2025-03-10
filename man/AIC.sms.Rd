% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{AIC.sms}
\alias{AIC.sms}
\alias{AIC}
\title{AIC}
\usage{
\method{AIC}{sms}(object, p = 2, n = Inf, ...)
}
\arguments{
\item{object}{Stock assesment object from smsR}

\item{p}{penalty on number of parameters, default = 2}

\item{n}{AICc sample size, default = Inf}

\item{...}{etc}
}
\value{
return the AIC value of the fitted assessment
}
\description{
Calculate AIC from smsR stock assessment object
}
\examples{

AIC(sas)
}
