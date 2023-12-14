prop_imp_quick <- function(f, CI_list, conf_vec=NULL, conf_final=0.95,...){

  K <- length(CI_list)
  if(is.null(conf_vec)) conf_vec <- rep(0.95, K)
  z_vec <- qnorm(1-(1-conf_vec)/2)
  SE_vec <- numeric(K)
  for(i in 1:K) SE_vec[i] <- diff(CI_list[[i]])/(2*z_vec[i])
  est_vec <- numeric(K)
  for(i in 1:K) est_vec[i] <- mean(CI_list[[i]])
  x <- array(est_vec,dim=c(K,1))
  deriv <- madness::numderiv(function(x){f(x,...)},x)
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
  return(list(est=f(est_vec,...),lower=f(as.numeric(CI_lower_new),...),upper=f(as.numeric(CI_upper_new),...)))
}

paf_levin_int <- function(x, thedim){
  prev <- c(1-sum(x[1:thedim]),x[1:thedim])
  lRR <- c(0,x[(thedim+1):(2*thedim)])
  (sum(exp(lRR)*prev)-1)/sum(exp(lRR)*prev)
}

paf_levin_cor_int <- function(x, thedim){
  prev <- c(1-sum(x[1:thedim]),x[1:thedim])
  lRR <- c(0,x[(thedim+1):(2*thedim)])
  lRRu <- c(0,x[(2*thedim+1):(3*thedim)])
  K <- thedim+1
  return(sum(prev[2:K]*exp(lRRu[2:K])*(exp(lRR[2:K])-1)/exp(lRR[2:K]))/(prev[1]+sum(prev[2:K]*exp(lRRu)[2:K])))
}

#' Implementation of Levin's formula for summary data
#'
#' @param prev A vector of estimated prevalence for each non-reference level of risk factor.  Can be left unspecified if conf_prev specified.
#' @param RR  A vector of estimated relative risk for each non-reference level of risk factor.  Can be left unspecified if conf_RR specified.
#' @param conf_prev If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for prevalence.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for prevalence of each non-reference level of risk factor.
#' @param conf_RR If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for relative risk.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for relative risk for each non-reference level of risk factor
#' @param digits integer.  The number of significant digits for rounding of PAF estimates and confidence intervals.  Default of 3.
#' @references Ferguson, J., Alvarez-Iglesias, A., Mulligan, M., Judge, C. and O’Donnell, M., 2024. Bias assessment and correction for Levin's population attributable fraction under confounding.  European Journal of Epidemiology, In press
#' @return If confidence intervals for prevalence and relative risk are not specified, the estimated PAF.  If confidence intervals for prevalence and relative risk are specified, confidence intervals for PAF are estimated using approximate propagation of imprecision.  Note that if confidence intervals are supplied as arguments, the algorithm makes assumptions that the point estimate of prevalence is the average of the specified confidence limits for prevalence, that the point estimate for relative risk is the geometric mean of the confidence limits for relative risk, and that the 3 estimators are independent.
#' @export
#' @examples
#' CI_p <- c(0.1,0.3)
#' CI_RR <- c(1.2, 2)
#' # calculation without confidence interval
#' paf_levin(prev=0.2,RR=exp(.5*log(1.2)+.5*log(2)))
#' # calculation with confidence interval
#' paf_levin(conf_prev=CI_p,conf_RR=CI_RR)
#' # add another level to risk factor
#' # with higher prevalence and RR
#' # this will increase the PAF
#' CI_p <- matrix(c(0.1,0.3,0.15, 0.25),nrow=2)
#' CI_RR <- matrix(c(1.2,2,1.5,3),nrow=2)
#' paf_levin(conf_prev=CI_p,conf_RR=CI_RR)
paf_levin <- function(prev=NULL, RR=NULL, conf_prev=NULL, conf_RR=NULL, digits=3){

  if(is.null(conf_prev)) prev <- c(1-sum(prev),prev)
  if(is.null(conf_RR)) RR <- c(1, RR)

  if(is.null(conf_prev) || is.null(conf_RR)) return(round((sum(RR*prev)-1)/sum(RR*prev),digits))

  if(!is.matrix(conf_prev)) conf_prev <- matrix(conf_prev,ncol=2)
  if(!is.matrix(conf_RR)) conf_RR <- matrix(conf_RR,ncol=2)

  log_conf_RR <- log(conf_RR)
  mydim <- nrow(conf_prev)
  thelist <- list()
  for(i in 1:mydim) thelist[[i]] <- conf_prev[i,]
  for(i in (mydim+1):(2*mydim)) thelist[[i]] <- log_conf_RR[i-mydim,]


   output <- prop_imp_quick(paf_levin_int, CI_list=thelist, thedim=mydim)
  return(paste(round(output$est, digits)," (", round(output$lower,digits),",",round(output$upper,digits),")",sep=""))


}

#' Implementation of Miettinen's formula for summary data
#'
#' @param prev A vector of estimated prevalence for each non-reference of risk factor.  Can be left unspecified if conf_prev specified.
#' @param RR  A vector of estimated causal relative risk for each non-reference level of risk factor.  Can be left unspecified if conf_RR specified.
#' @param RRu  A vector of estimated unadjusted relative risk for each non-reference level of the risk factor.  Can be left unspecified if conf_RRu specified.
#' @param conf_prev If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for prevalence.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for prevalence of each non-refernece level.
#' @param conf_RR If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for the causal relative risk.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for causal relative risk fror each non-reference level of risk factor.
#' @param conf_RRu If risk factor has 2 levels, a numeric vector of length 2 giving confidence limits for the unadjusted relative risk.  If risk factor has K>2 levels, a K-1 x 2 matrix giving confidence intervals for unadjusted relative risk for each non-reference level of risk factor.
#' @param digits integer.  The number of significant digits for rounding of PAF estimates and confidence intervals.  Default of 3
#' @references Ferguson, J., Alvarez-Iglesias, A., Mulligan, M., Judge, C. and O’Donnell, M., 2024. Bias assessment and correction for Levin's population attributable fraction under confounding.  European Journal of Epidemiology, In press
#' @return If confidence intervals for prevalence, adjusted and unadjusted relative risk are not specified, the estimated PAF.  If confidence intervals are specified, confidence intervals for PAF are also estimated using approximate propagation of imprecision.  Note that if confidence intervals are supplied as arguments, the algorithm makes assumptions that the point estimate of prevalence is the average of the specified confidence limits for prevalence, the point estimates for adjusted/unadjusted relative risk are the geometric means of the specified confidence limits for relative risk, and that the 3 estimators are independent.
#' @export
#' @examples
#' CI_p <- c(0.1,0.3)
#' CI_RR <- c(1.2, 2)
#' CI_RRu <- c(1.5, 2.5)
#' # example without confidence interval
#' paf_miettinen(prev=0.2,RR=exp(.5*log(1.2)+.5*log(2)), RRu=exp(.5*log(1.5)+.5*log(2.5)))
#' #' # example with confidence interval
#' paf_miettinen(conf_prev=CI_p,conf_RR=CI_RR, conf_RRu=CI_RRu)
#' # risk factor with more than two non-reference levels
#' # confidence intervals for non-reference levels
#' # of risk factor should be a (K-1) x 2 matrix
#' CI_p <- matrix(c(0.1,0.3,0.15, 0.25),nrow=2)
#' CI_RR <- matrix(c(1.2,2,1.5,3),nrow=2)
#' CI_RRu <- matrix(c(1.5,2.5,2,3.5),nrow=2)
#' paf_miettinen(conf_prev=CI_p,conf_RR=CI_RR, conf_RRu=CI_RRu)
paf_miettinen  <- function(prev=NULL, RR=NULL, RRu=NULL, conf_prev=NULL, conf_RR=NULL, conf_RRu=NULL,digits=3){

  if(!is.null(prev)) prev <- c(1-sum(prev),prev)
  if(!is.null(RR)) lRR <- c(0, log(RR))
  if(!is.null(RRu)) lRRu <- c(0, log(RRu))
  K <- length(RR)+1

  if(is.null(conf_prev) || is.null(conf_RR) || is.null(conf_RRu)) return(round(sum(prev[2:K]*exp(lRRu[2:K])*(exp(lRR[2:K])-1)/exp(lRR[2:K]))/(prev[1]+sum(prev[2:K]*exp(lRRu)[2:K])),digits))

  if(!is.matrix(conf_prev)) conf_prev <- matrix(conf_prev,ncol=2)
  if(!is.matrix(conf_RR)) conf_RR <- matrix(conf_RR,ncol=2)
  if(!is.matrix(conf_RRu)) conf_RRu <- matrix(conf_RRu,ncol=2)


  log_conf_RR <- log(conf_RR)
  log_conf_RRu <- log(conf_RRu)
  mydim <- nrow(conf_prev)
  thelist <- list()
  for(i in 1:mydim) thelist[[i]] <- conf_prev[i,]
  for(i in (mydim+1):(2*mydim)) thelist[[i]] <- log_conf_RR[i-mydim,]
  for(i in (2*mydim+1):(3*mydim)) thelist[[i]] <- log_conf_RRu[i-2*mydim,]

  output <- prop_imp_quick(paf_levin_cor_int, CI_list=thelist, thedim=mydim)
  return(paste(round(output$est, digits)," (", round(output$lower,digits),",",round(output$upper,digits),")",sep=""))
}

CI_p <- c(0.1,0.3)
CI_RR <- c(1.2, 2)
CI_RRu <- c(1.5, 2.5)
