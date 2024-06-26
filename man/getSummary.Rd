% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-output.R
\name{getSummary}
\alias{getSummary}
\title{Get summary of derived variables from fitted smsR object}
\usage{
getSummary(df.tmb, sas)
}
\arguments{
\item{df.tmb}{list of input parameters for smsR object}

\item{sas}{fitted smsR model}
}
\value{
returns a table of SSB, R, Yield and Fbar with uncertainty
}
\description{
Get summary of derived variables from fitted smsR object
}
\examples{

dat <- getSummary(df.tmb, sas)

}
