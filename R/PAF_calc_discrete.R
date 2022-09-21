#' Calculation of attributable fraction using a categorized risk factor
#'
#' @param model Either a clogit, glm or coxph fitted regression object.  Non-linear effects can be specified in these models if necessary via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param riskfactor The name of the risk factor of interest in the dataset.  The risk factor values can be 0/1 numeric, categorical or factor valued.
#' @param refval The reference value for the risk factor.  If a risk factor is 0/1 numeric, 0 is assumed as the default value, otherwise refval must be specified.
#' @param data A dataframe containing variables used for fitting the model
#' @param calculation_method A character either 'B' (Bruzzi) or 'D' (Direct method).  For case control data, the method described in Bruzzi 1985 is recommended.  Bruzzi's method estimates PAF from relative risks and prevalence of exposure to the risk factor.  The Direct method estimates PAF by summing estimated probabilities of disease in the absence of exposure on the individual level
#' @param prev The estimated prevalence of disease (A number between 0 and 1).  This only needs to be specified if the data source is from a case control study, and the direct method is used
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with point estimate (default FALSE)
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param t_vector Numeric.  A vector of times at which to calculate PAF (only specified if model is coxph)
#' @param ci_level Numeric.  A number between 0 and 1 specifying the confidence level
#' @param ci_type Character.  A vector specifying the types of confidence interval desired, as available in the 'Boot' package. The default is c('norm'), which calculates a symmetric confidence interval: (Est-Bias +- 1.96*SE), with the standard error calculated via Bootstrap.  Other choices are 'basic', 'perc' and 'bca'.  Increasing the number of Bootstrap repetitions is recommended for the 'basic', 'perc' and 'bca' methods.
#' @references Bruzzi, P., Green, S.B., Byar, D.P., Brinton, L.A. and Schairer, C., 1985. Estimating the population attributable risk for multiple risk factors using case-control data. American journal of epidemiology, 122(5), pp.904-914
#' @return An estimated PAF if ci=FALSE, or for survival data a vector of estimated PAF corresponding to event times in the data.  If ci=TRUE, a vector with elements corresponding to the raw estimate, estimated bias, bias corrected estimate and lower and upper elements of any confidence procedures requested.  If ci=TRUE, and a coxph model is fit, a matrix will be returned, with rows corresponding to differing times at which the PAF might be calculated.
#' @export
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' data(stroke_reduced)
#' model_exercise <- glm(formula = case ~ region * ns(age, df = 5) +
#'  sex * ns(age, df = 5) + education + exercise + ns(diet, df = 3) +
#'  smoking + alcohol + stress, family = "binomial", data = stroke_reduced)
#' # calculate discrete PAF using Bruzzi method
#' PAF_calc_discrete(model_exercise, "exercise", refval=0,
#' data=stroke_reduced, calculation_method="B",ci=FALSE)
#' \donttest{
#' # calculate discrete PAF using Direct method
#' # PAF_calc_discrete(model_exercise, "exercise", refval=0,
#' # data=stroke_reduced, calculation_method="D", prev=0.0035)
#' # Use bootstrap resampling to calculate a confidence interval
#' PAF_calc_discrete(model_exercise, "exercise", refval=0,
#' data=stroke_reduced, calculation_method="D", prev=0.005, ci=TRUE, boot_rep=10)
#' ### use the Bruzzi method derived by Bruzzi, 1985, instead
#' PAF_calc_discrete(model_exercise, "exercise", refval=0, data=stroke_reduced,
#'  calculation_method="B", ci=TRUE, boot_rep=10)
#' # examples of clogit and coxph regressions
#'
#' model_high_blood_pressure_clogit <- clogit(formula = case ~ age +
#' education +exercise + ns(diet, df = 3) + smoking + alcohol + stress +
#'  ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) + high_blood_pressure +
#'  strata(strata),data=stroke_reduced)
#' PAF_calc_discrete(model_high_blood_pressure_clogit, "high_blood_pressure",
#' refval=0, data=stroke_reduced, calculation_method="B",ci=TRUE, boot_rep=100,
#'  ci_type=c('norm'))
#'
#' model_high_blood_pressure_coxph <- coxph(formula = Surv(time,event) ~
#'  ns(age,df=5) + education +exercise + ns(diet, df = 3) + smoking + alcohol +
#'   stress + ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) +
#'   high_blood_pressure, data = stroke_reduced)
#' PAF_calc_discrete(model_high_blood_pressure_coxph, "high_blood_pressure",
#' refval=0, data=stroke_reduced, calculation_method="D", ci=TRUE,
#' boot_rep=10, ci_type=c('norm'),t_vector=c(1,2,3,4,5,6,7,8,9))
#' }
PAF_calc_discrete <- function(model, riskfactor, refval, data, calculation_method="B", prev=NULL,ci=FALSE,boot_rep=100, t_vector=NULL, ci_level=.95, ci_type=c("norm")){

  if(!is.data.frame(data)){
      stop(
        "Data must be a dataframe object")
  }

  # remove data not used to fit model
  data <- data[row.names(data) %in% row.names(model.frame(model)),]

  new_data <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data=data)


  which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

  N <- nrow(data)
  riskfactor_vals <- data[,which_col]
  if(is.integer(riskfactor_vals)) data[,which_col] <- as.numeric(data[,which_col])

  # call the impact fraction function using the predicted data frame
  impact_fraction(model=model, data=data, new_data=new_data,calculation_method=calculation_method, prev=prev,ci=ci,boot_rep=boot_rep,t_vector=t_vector, ci_level= ci_level, ci_type=ci_type)
}





#' Create a data frame for predictions (when risk factor is discrete).
#'
#' @param riskfactor The name of the risk factor of interest in the dataset
#' @param refval The reference value for the risk factor
#' @param data A dataframe containing variables used to fit the model
#'
#' @return A data frame where the categorical variable is set to it's reference level
predict_df_discrete <- function(riskfactor, refval, data){

  if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

    stop("Riskfactor not in dataset.  Check spelling")

  }

  which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

  N <- nrow(data)
  riskfactor_vals <- data[,which_col]

  if(is.numeric(riskfactor_vals)){
    if(is.na(refval)) refval <- 0
    if(!all(riskfactor_vals %in% c(0,1)) || refval !=0){
      stop("Numeric risk factors must be 0/1, with the reference set to 0")

    }

  }

  if(is.numeric((data[,which_col]))) data[,which_col] <- rep(refval,N)
  if(is.character(data[,which_col])) data[,which_col] <- as.character(rep(refval,N))
  if(is.factor(data[,which_col])) data[,which_col] <- factor(rep(refval,N),levels = levels(data[,which_col]))

  return(data)

}





