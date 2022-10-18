#' Calculation of attributable fractions with a continuous exposure
#'
#' @param model Either a clogit, glm or coxph R model object.  Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param riskfactor_vec A character vector of names for continuous exposures/riskfactors that are predictors the model.
#' @param q_vec A vector of 'risk quantiles' for the continuous exposure.  q_vec=c(0.01) (the default) calculates an estimate of the PAF that is in some way analogous to eliminating a categorical risk factor.  Other values in q_vec correspond to interventions on the continuous risk factor that results in a risk level for all individuals thresholded above by the corresponding quantile of pre-intervention population risk.  For survival regressions only single element values of q_vec are allowed
#' @param data A dataframe containing variables used for fitting the model
#' @param calculation_method A character either 'B' (Bruzzi) or 'D' (Direct method).  For case control data, the method described in Bruzzi 1985 is recommended.  Bruzzi's method estimates PAF from relative risks and prevalence of exposure to the risk factor.  The Direct method estimates PAF via summing estimated probabilities of disease in the absence of exposure over differing individuals.
#' @param prev The estimated prevalence of disease (A number between 0 and 1).  This only needs to be specified if the data source is from a case control study, and the direct calculation method is used
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with point estimate (default FALSE)
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param t_vector Numeric.  A vector of times at which to calculate PAF (only specified if model is coxph)
#' @param ci_level Numeric.  A number between 0 and 1 specifying the confidence level
#' @param ci_type Character.  A vector specifying the types of confidence interval desired, as available in the 'Boot' package. The default is c('norm'), which calculates a symmetric confidence interval: (Est-Bias +- 1.96*SE), with the standard error calculated via Bootstrap.  Other choices are 'basic', 'perc' and 'bca'.  Increasing the number of Bootstrap repetitions is recommended for the 'basic', 'perc' and 'bca' methods.
#' @param S Integer (default 1).  Only relevant to change if there is an interaction between the continuous exposure and other variables in the model.  In this case, marginal comparisons of disease risk at differing levels of the exposure need to be averaged over many individuals.  S is the number of individuals used in this averaging.  May be slow for large S
#' @param weight_vec An optional vector of inverse sampling weights for survey data (note that variance will not be calculated correctly if sampling isn't independent).  Note that this vector will be ignored if prev is specified, and the weights will be calibrated so that the weighted sample prevalence of disease equals prev.
#' @return A PAF_q object.  When ci=FALSE, this will essentially be a vector of estimated PAF corresponding to the quantiles specified in q_vec.  If ci=TRUE, a data frame with columns corresponding to the raw estimate, estimated bias, bias corrected estimate and lower and upper elements of any confidence procedures requested, and rows corresponding to the quantiles in q_vec.
#' @export
#'
#' @references Ferguson, J., Maturo, F., Yusuf, S. and O’Donnell, M., 2020. Population attributable fractions for continuously distributed exposures. Epidemiologic Methods, 9(1).
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' # Example with logistic regression.  PAF_q (as in Ferguson, 2020)
#' # estimated at q=0.01, 0.1, 0.3, 0.5, 0.7. 0.9.  PAF_0.01 is roughly
#' # analogous to 'eliminating' a discrete risk factor, but its estimation
#' # may be unstable for some exposures, and the corresponding intervention
#' # may be impractical.  Comparing PAF_q for q >= 0.1 over different risk factors
#' # may lead to more sensible comparisons of disease burden.
#' # Either method (direct, D, or Bruzzi )
#' # reduce dataset to improve run time (not recommended on real data!)
#' stroke_small <- stroke_reduced[sample(1:nrow(stroke_reduced),1000),]
#' model_continuous <- glm(formula = case ~ region * ns(age, df = 5) +
#' sex * ns(age, df = 5) + education +exercise + ns(diet, df = 3) +
#' alcohol + stress + ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) +
#'  high_blood_pressure, family = "binomial", data = stroke_small)
#' out <- PAF_calc_continuous(model_continuous,riskfactor_vec=
#' c("diet","lipids","waist_hip_ratio"),q_vec=c(0.1,0.5,0.9),
#' ci=FALSE,calculation_method="D",data=stroke_small, prev=0.0035)
#' print(out)
#' plot(out)
#' \donttest{
#'  # with confidence intervals (via bootstrap) on full dataset.  Slower.
#' model_continuous_clogit <- clogit(formula = case ~ region * ns(age, df = 5) +
#' sex * ns(age, df = 5) + education +exercise + ns(diet, df = 3)  +
#' alcohol + stress + ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) +
#'  high_blood_pressure + strata(strata), data = stroke_reduced)
#' out <- PAF_calc_continuous(model_continuous_clogit,riskfactor_vec=c("diet",
#' "lipids","waist_hip_ratio"),q_vec=c(0.01, 0.1,0.3,0.5,0.7,0.9),
#' ci=TRUE,calculation_method="B",data=stroke_reduced, prev=0.01)
#' print(out)
#' plot(out)
#' }
PAF_calc_continuous <- function(model, riskfactor_vec, q_vec=c(0.01), data, calculation_method="B", prev=NULL,ci=FALSE,boot_rep=10, t_vector=NULL, ci_level=.95, ci_type=c("norm"), S=1, weight_vec=NULL){

  if(is.null(weight_vec)) weight_vec <- rep(1,nrow(data))
  if(!is.data.frame(data)){
    stop(
      "Data must be a dataframe object")
  }
  data <- as.data.frame(data) # in case a tibble which passes the above test

  # remove data not used to fit model
  data <- data[row.names(model.frame(model)) %in% row.names(data),]

  model_type <- class(model)[1]
  if(model_type=="glm"){
      if (!as.character(model$family[1])=="binomial" & ! as.character(model$family[2]) %in% c("logit","log")) {
      stop(
        "The family must be binomial and link must be either log or logistic"
      )
    }
  }

  N <- nrow(data)

  if(model_type=="coxph" && length(t_vector)>1){

    stop(
      "Function can calculate vector of PAF_q for continuous risk factors, only at a single time point"
    )

  }


   if(!ci) {
      res <- impact_fraction_qvec(data, ind=(1:N), model=model, model_type=model_type, riskfactor_vec=riskfactor_vec,  q_vec=q_vec,calculation_method=calculation_method,S=S,prev=prev,t_vector=t_vector, weight_vec=weight_vec)
      q_vec_obj <- structure(data.frame("riskfactor"=rep(riskfactor_vec,times=rep(length(q_vec),length(riskfactor_vec))),"q_val"=rep(q_vec,length(riskfactor_vec)),paf_q=res),class="PAF_q")

         return(q_vec_obj)

   }
    if(ci) {
      nc <- options()$boot.ncpus
      cl <- parallel::makeCluster(nc)
      parallel::clusterExport(cl, c("coxph","clogit","strata","Surv"))
      parallel::clusterExport(cl, c("ns"))
      parallel::clusterExport(cl, c("predict_df_continuous","risk_quantiles","impact_fraction","if_bruzzi","if_direct"))
      res <- boot::boot(data, statistic=impact_fraction_qvec,model=model,  model_type=model_type, riskfactor_vec=riskfactor_vec, q_vec=q_vec,calculation_method=calculation_method,S=S, prev=prev,t_vector=t_vector, R=boot_rep,cl=cl, weight_vec=weight_vec)
      parallel::stopCluster(cl)
      q_vec_obj <- structure(cbind("riskfactor"=rep(riskfactor_vec,times=rep(length(q_vec),length(riskfactor_vec))),"q_val"=rep(q_vec,length(riskfactor_vec)),extract_ci(res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type,continuous=TRUE)),class="PAF_q")
            return(q_vec_obj)
    }

}


##  write functions that that risk_factor vectors and q_vec and risk_q lists and give corresponding impact fractions
impact_fraction_qvec <- function(data, ind, model, model_type, riskfactor_vec,  q_vec, calculation_method, S=1, prev=NULL,t_vector=NULL, weight_vec=NULL){

  N <- nrow(data)

  resamples <- (1:N)
  if(model_type == "clogit"){

    if(!all(ind==(1:N))){

      model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
      model_text <- paste0(model_text[2],model_text[1],model_text[3])
      strataname <- gsub(".*strata\\((.*)\\).*",replacement="\\1",x=model_text,perl=TRUE)

      # find strata variable
      strataids <- data[,colnames(data)==strataname]
      validids <- names(table(strataids))[table(strataids)==2]
      possibleresamples <- (1:nrow(data))[strataids %in% validids]
      ## order possible resamples indexes according to valid ids
      possibleresamples <- possibleresamples[order(strataids[strataids %in% validids])]
      totake <- sample(1:(length(possibleresamples)/2),length(possibleresamples)/2,replace=TRUE)
      resamples <- c(possibleresamples[2*totake],possibleresamples[2*totake-1])
      data <- data[resamples,]
      # avoid duplication of strata names
      data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))


      #refit model
      model_text <- paste0("survival::clogit(",model_text,",data=data)")
      thesplit <- ""
      while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
        model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
        stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
        model_text <- stuff[[1]][1]
        thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
      }
      model_text <- paste0(model_text,thesplit)

      model <- eval(parse(text=model_text))

    }
  }
  if(model_type == "coxph"){

    if(!all(ind==(1:N))){

      data <- data[ind, ]
      if(!is.null(weight_vec)) weight_vec <- weight_vec[ind]
      model <- update(model,data=data)
    }
  }
  if(model_type == "glm"){

    if(!all(ind==(1:N))){

      data <- data[ind, ]
      if(!is.null(weight_vec)) weight_vec <- weight_vec[ind]
      model <- update(model,data=data)
    }
  }
    risk_q <- list()
    for(i in 1:length(riskfactor_vec)){
      risk_q[[i]] <- risk_quantiles(riskfactor=riskfactor_vec[i], data=data, model=model, S=S, q=seq(from=0.01,to=0.99,by=0.01))
    }

   return_vec <- c()
    for(i in 1:length(riskfactor_vec)){
        for(j in 1:length(q_vec)){

          new_data <- predict_df_continuous(riskfactor=riskfactor_vec[i], q_val=q_vec[j],risk_q=risk_q[[i]], data=data)
          if(i==1 & j==1) return_vec <- impact_fraction(model=model, data=data, new_data=new_data, calculation_method=calculation_method,prev=prev,ci=FALSE,t_vector=t_vector, weight_vec = weight_vec)
          if(i!=1 || j!=1) return_vec <- c(return_vec,impact_fraction(model=model, data=data, new_data=new_data, calculation_method=calculation_method,prev=prev,ci=FALSE,t_vector=t_vector, weight_vec=weight_vec))
        }
    }
  return(return_vec)
}



#' Internal:  Create a data frame for predictions (when risk factor is continuous).
#'
#' @param riskfactor The name of the risk factor of interest in the dataset
#' @param q_val The risk quantile to match to
#' @param risk_q Estimated risk quantiles
#' @param data A dataframe containing variables used to fit the model
#' @return A data frame where the distribution continuous risk factor so at an individual level, risk is at the q_val-quantile or below
#' @export
predict_df_continuous <- function(riskfactor, q_val,risk_q, data){

  if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

    stop("Riskfactor not in dataset.  Check spelling")

  }

  which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

  N <- nrow(data)
  riskfactor_vals <- data[,which_col]

  if(!is.numeric(riskfactor_vals)){

    stop("Risk factor values must be numeric")

  }
  index_vec <- sapply(riskfactor_vals,function(x){which.min((x-risk_q)^2)})
  est_quantiles <- as.numeric(names(risk_q)[index_vec])
  intervention <- risk_q[which.min((as.numeric(names(risk_q))-q_val)^2)]
  riskfactor_vals[est_quantiles >= q_val] <- intervention
  data[,which_col] <- riskfactor_vals

  return(data)

}



#' Return the vector of risk quantiles for a continuous risk factor.
#'
#' @param riskfactor The name of the risk factor of interest in the dataset.  A character vector
#' @param data A dataframe containing variables used to fit the model
#' @param model The fitted model
#' @param q The desired risk quantiles
#' @param S The number of randomly selected individuals for which risk is measured (defaults to 1).  Let to perhaps 100 if risk factor involved in interactions in model
#'
#' @return A named vector of size S giving the risk factor quantiles
#' @export
risk_quantiles <- function(riskfactor, data, model, S=1, q=seq(from=0.01,to=0.99,by=0.01)){

  data <- data[row.names(model.frame(model)) %in% row.names(data),]

  if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

    stop("Riskfactor not in dataset.  Check spelling")

  }

  which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

  N <- nrow(data)
  riskfactor_vals <- data[,which_col]
  if(!is.numeric(riskfactor_vals)){

    stop("Risk factor variable must be numeric")

  }

  subset <- sample(1:N,size=S,replace=TRUE)
  newdata <- data
  predict_vals_mat <- matrix(nrow=N,ncol=S)
  for(i in 1:S){

    newdata[,1:ncol(data)] <- matrix(rep(as.character(data[subset[i],]),N),nrow=N,byrow=TRUE)

    for(j in 1:ncol(newdata)){

      if(is.factor(data[,j])) newdata[,j] <- factor(newdata[,j],levels=levels(data[,j]))
      if(is.numeric(data[,j])) newdata[,j] <- as.numeric(newdata[,j])
    }
    newdata[,which_col] <- riskfactor_vals
    predict_vals_mat[,i] <- predict(model,newdata=newdata)

  }
  predict_vals_summary <- apply(predict_vals_mat,1,mean,na.rm=TRUE)
  return_vec <- riskfactor_vals[order(predict_vals_summary)][ceiling(N*q)]
  names(return_vec) <- q
  return(return_vec)
}

#' Print out PAF_q for differing risk factors
#'
#' @param x A PAF_q object.  This is a dataframe that is created by running the function PAF_calc_continuous. The final 3 columns of the data frame are assumed to be (in order), PAF and lower and upper confidence bounds.
#' @param ...  Other arguments to be passed to print
#' @return No return value, prints the PAF_q object to the console.
#' @export
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' model_continuous <- glm(formula = case ~ region * ns(age, df = 5) +
#' sex * ns(age, df = 5) + education +exercise + ns(diet, df = 3) +
#'  alcohol + stress + ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) +
#' high_blood_pressure, family = "binomial", data = stroke_reduced)
#' out <- PAF_calc_continuous(model_continuous,
#' riskfactor_vec=c("diet","lipids","waist_hip_ratio"),
#' q_vec=c(0.01, 0.1,0.3,0.5,0.7,0.9),ci=FALSE,calculation_method="B",
#' data=stroke_reduced)
#' print(out)
print.PAF_q <- function(x,...){

  data_frame <- structure(as.list(x),class="data.frame", row.names=attr(x,"row.names"))
  print(data_frame)

}


