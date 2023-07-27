prop_imp_quick <- function(f, CI_list, conf_vec=NULL, conf_final=0.95,...){

  K <- length(CI_list)
  if(is.null(conf_vec)) conf_vec <- rep(0.95, K)
  z_vec <- qnorm(1-(1-conf_vec)/2)
  SE_vec <- numeric(K)
  for(i in 1:K) SE_vec[i] <- diff(CI_list[[i]])/(2*z_vec[i])
  est_vec <- numeric(K)
  for(i in 1:K) est_vec[i] <- mean(CI_list[[i]])
  x <- array(est_vec,dim=c(K,1))
  deriv <- madness::numderiv(f,x)
  abs_deriv <- abs(deriv) # we'll assume function is parameterised so that it is increasing in all arguments
  q_vec <- abs_deriv*SE_vec
  z <- qnorm(1-(1-conf_final)/2)
  z_star <- (z/sqrt(sum(q_vec^2)))*q_vec
  CI_upper_new <- est_vec + z_star*SE_vec
  CI_lower_new <- est_vec - z_star*SE_vec
  if(any(deriv<0)){
    CI_upper_new[deriv<0] <- est_vec[deriv<0] - z_star[deriv<0]*SE_vec[deriv<0]
    CI_lower_new[deriv<0] <- est_vec[deriv<0] + z_star[deriv<0]*SE_vec[deriv<0]
  }
  return(list(est=f(est_vec),lower=f(as.numeric(CI_lower_new)),upper=f(as.numeric(CI_upper_new))))
}

paf_levin_int <- function(x){
  prev <- x[1]
  lRR <- x[2]
  prev*(exp(lRR)-1)/(1+prev*(exp(lRR)-1))
}

paf_levin_cor_int <- function(x){
  prev <- x[1]
  lRR <- x[2]
  lRRu <- x[3]
  return(prev*exp(lRRu)/(1+prev*(exp(lRRu)-1))*((exp(lRR)-1)/exp(lRR)))
}

#' Levin's formula based on relative risk and prevalence
#'
#' @param prev Estimated prevalence.  Can be left unspecified if conf_prev specified.
#' @param RR Estimated relative risk.  Can be left unspecified if conf_RR specified.
#' @param conf_prev A numeric vector of length 2 giving confidence limits for prevalence.
#' @param conf_RR A numeric vector of length 2 giving the confidence limits for relative risk.#'
#' @param digits integer.  The number of significant digits for rounding of PAF estimates and confidence intervals.  Default of 3.
#' @return If confidence intervals for prevalence and relative risk are not specified, the estimated PAF.  If confidence intervals for prevalence and relative risk are specified, confidence intervals for PAF are estimated using approximate propagation of imprecision.  Note that if confidence intervals are supplied as arguments, the algorithm makes assumptions that the point estimate of prevalence is the average of the specified confidence limits for prevalence, the point estimate for relative risk is the geometric mean of the confidence limits for relative risk, and that the 3 estimators are independent.
#' @export
#' @examples
#' CI_p <- c(0.1,0.3)
#' CI_RR <- c(1.2, 2)
#' # calculation without confidence interval
#' paf_levin(prev=0.2,RR=exp(.5*log(1.2)+.5*log(2)))
#' # calculation with confidence interval
#' paf_levin(conf_prev=CI_p,conf_RR=CI_RR)
paf_levin <- function(prev=NULL, RR=NULL, conf_prev=NULL, conf_RR=NULL, digits=3){

  if(is.null(conf_prev) || is.null(conf_RR)) return(round(prev*(RR-1)/(1+prev*(RR-1)),digits))

  if(!is.numeric(conf_prev) || !length(conf_prev==2) || !is.numeric(conf_RR) || !length(conf_RR==2)) return("Error: confidence intervals should be numeric vectors of length 2")

  CI_RR <- log(conf_RR)
   output <- prop_imp_quick(paf_levin_int, CI_list=list(conf_prev,CI_RR))
  return(paste(round(output$est, digits)," (", round(output$lower,digits),",",round(output$upper,digits),")",sep=""))


}

#' Miettinen's formula based on adjusted relative risk, unadjusted relative risk and prevalence
#'
#' @param prev Estimated prevalence.  Can be left unspecified if conf_prev specified.
#' @param RR Estimated adjusted relative risk.  Can be left unspecified if conf_RR specified.
#' @param RR_u Estimated unadjusted relative risk.  Can be left unspecified if conf_RRu specified.
#' @param conf_prev A numeric vector of length 2 giving confidence limits for prevalence.
#' @param conf_RR A numeric vector of length 2 giving the confidence limits for the adjusted relative risk.
#' @param conf_RRu A numeric vector of length 2 giving the confidence limits for the unadjusted relative risk.
#' @param digits integer.  The number of significant digits for rounding of PAF estimates and confidence intervals.  Default of 3.
#' @return If confidence intervals for prevalence and adjusted and unadjusted relative risk are not specified, the estimated PAF.  If confidence intervals are specified, confidence intervals for PAF are also estimated using approximate propagation of imprecision.  Note that if confidence intervals are supplied as arguments, the algorithm makes assumptions that the point estimate of prevalence is the average of the specified confidence limits for prevalence, the point estimates for adjusted/unadjusted relative risk are the geometric means of the specified confidence limits for relative risk, and that the 3 estimators are independent.
#' @export
#' @examples
#' CI_p <- c(0.1,0.3)
#' CI_RR <- c(1.2, 2)
#' CI_RRu <- c(1.5, 2.5)
#' # example without confidence interval
#' paf_miettinen(prev=0.2,RR=exp(.5*log(1.2)+.5*log(2)), RR_u=exp(.5*log(1.5)+.5*log(2.5)))
#' #' # example with confidence interval
#' paf_miettinen(conf_prev=CI_p,conf_RR=CI_RR, conf_RRu=CI_RRu)
paf_miettinen  <- function(prev=NULL, RR=NULL, RR_u=NULL, conf_prev=NULL, conf_RR=NULL, conf_RRu=NULL,digits=3){


  if(is.null(conf_prev) || is.null(conf_RR) || is.null(conf_RRu)) return(round(prev*RR_u/(1+prev*(RR_u-1))*(RR-1)/RR,digits))

  if(!is.numeric(conf_prev) || !length(conf_prev==2) || !is.numeric(conf_RRu) || !length(conf_RRu==2) || !is.numeric(conf_RR) || !length(conf_RR==2)) return("Error: confidence intervals should be numeric vectors of length 2")

  CI_lRR <- log(conf_RR)
  CI_lRRu <- log(conf_RRu)
  output <- prop_imp_quick(paf_levin_cor_int, CI_list=list(conf_prev,CI_lRR,CI_lRRu))
  return(paste(round(output$est, digits)," (", round(output$lower,digits),",",round(output$upper,digits),")",sep=""))
}

CI_p <- c(0.1,0.3)
CI_RR <- c(1.2, 2)
CI_RRu <- c(1.5, 2.5)
