% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Levin.R
\name{paf_miettinen}
\alias{paf_miettinen}
\title{Implementation of Miettinen's formula for summary data}
\usage{
paf_miettinen(
  prev = NULL,
  RR = NULL,
  RRu = NULL,
  conf_prev = NULL,
  conf_RR = NULL,
  conf_RRu = NULL,
  digits = 3
)
}
\arguments{
\item{prev}{A vector of estimated prevalence for each non-reference of risk factor.  Can be left unspecified if conf_prev specified.}

\item{RR}{A vector of estimated causal relative risk for each non-reference level of risk factor.  Can be left unspecified if conf_RR specified.}

\item{RRu}{A vector of estimated unadjusted relative risk for each non-reference level of the risk factor.  Can be left unspecified if conf_RRu specified.}

\item{conf_prev}{If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for prevalence.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for prevalence of each non-refernece level.}

\item{conf_RR}{If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for the causal relative risk.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for causal relative risk fror each non-reference level of risk factor.}

\item{conf_RRu}{If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for the unadjusted relative risk.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for unadjusted relative risk for each non-reference level of risk factor.}

\item{digits}{integer.  The number of significant digits for rounding of PAF estimates and confidence intervals.  Default of 3}
}
\value{
If confidence intervals for prevalence, adjusted and unadjusted relative risk are not specified, the estimated PAF.  If confidence intervals are specified, confidence intervals for PAF are also estimated using approximate propagation of imprecision.  Note that if confidence intervals are supplied as arguments, the algorithm makes assumptions that the point estimate of prevalence is the average of the specified confidence limits for prevalence, the point estimates for adjusted/unadjusted relative risk are the geometric means of the specified confidence limits for relative risk, and that the 3 estimators are independent.
}
\description{
Implementation of Miettinen's formula for summary data
}
\examples{
CI_p <- c(0.1,0.3)
CI_RR <- c(1.2, 2)
CI_RRu <- c(1.5, 2.5)
# example without confidence interval
paf_miettinen(prev=0.2,RR=exp(.5*log(1.2)+.5*log(2)), RRu=exp(.5*log(1.5)+.5*log(2.5)))
#' # example with confidence interval
paf_miettinen(conf_prev=CI_p,conf_RR=CI_RR, conf_RRu=CI_RRu)
# risk factor with more than two non-reference levels
# confidence intervals for non-reference levels
# of risk factor should be a (K-1) x 2 matrix
CI_p <- matrix(c(0.1,0.3,0.15, 0.25),nrow=2)
CI_RR <- matrix(c(1.2,2,1.5,3),nrow=2)
CI_RRu <- matrix(c(1.5,2.5,2,3.5),nrow=2)
paf_miettinen(conf_prev=CI_p,conf_RR=CI_RR, conf_RRu=CI_RRu)
}
\references{
Ferguson, J., Alvarez-Iglesias, A., Mulligan, M., Judge, C. and O’Donnell, M., 2024. Bias assessment and correction for Levin's population attributable fraction under confounding.  European Journal of Epidemiology, In press
}
