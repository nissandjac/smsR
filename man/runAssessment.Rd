% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runAssessment.R
\name{runAssessment}
\alias{runAssessment}
\title{run the assessment with the data and the estimated parameters}
\usage{
runAssessment(
  df.tmb,
  parms,
  lwr = list(NA),
  upr = list(NA),
  mps = NULL,
  silent = TRUE,
  dll = "smsR"
)
}
\arguments{
\item{df.tmb}{list of input data}

\item{parms}{parameters}

\item{lwr}{optional list of user evaluated lower bounds}

\item{upr}{optional list of user evaluated lower bounds}

\item{mps}{optional list of parameteres to be mapped}

\item{silent}{silent or verbose}

\item{dll}{optional file in case of custom stuff}
}
\description{
run the assessment with the data and the estimated parameters
}
