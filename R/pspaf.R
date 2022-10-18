#' Estimate pathway specific population attributable fractions
#'
#' @param response_model A R model object for a binary outcome that involves a risk factor, confounders and mediators of the risk factor outcome relationship.  Note that a weighted model should be used for case control data.  Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param mediator_models A list of fitted  models describing the risk factor/mediator relationship (the predictors in the model will be the risk factor and any confounders)  Note a weighted model should be fit when data arise from a case control study.  Models can be specified for linear responses (lm), binary responses (glm) and ordinal factors (through polr).  Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param riskfactor character.  Represents the name of the risk factor
#' @param refval For factor valued risk factors, the reference level of the risk factor.  If the risk factor is numeric, the reference level is assumed to be 0.
#' @param data dataframe. A dataframe (with no missing values) containing the data used to fit the mediator and response models.  You can run data_clean to the input dataset if the data has missing values as a pre-processing step
#' @param prev numeric.  A value between 0 and 1 specifying the prevalence of disease: only relevant to set if data arises from a case control study.
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param ci logical.  If TRUE a confidence interval is calculated using Bootstrap
#' @param ci_level Numeric.  Default 0.95. A number between 0 and 1 specifying the confidence level (only necessary to specify when ci=TRUE)
#' @param ci_type Character.  Defalt norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @param weight_vec An optional vector of inverse sampling weights for survey data (note that variance will not be calculated correctly if sampling isn't independent).  Note that this will be ignored if prev is specified and calculation_method="D", in which case the weights will be constructed so the empirical re-weighted prevalence of disease is equal to prev
#' @return A numeric vector (if ci=FALSE), or data frame (if CI=TRUE) containing estimated PS-PAF for each mediator referred to in mediator_models, together with estimated direct PS-PAF and possibly confidence intervals.
#' @export
#'
#' @references Pathway specific Population attributable fractions.  O’Connell, M.M. and Ferguson, J.P., 2022. IEA. International Journal of Epidemiology, 1, p.13.  Accessible at: https://academic.oup.com/ije/advance-article/doi/10.1093/ije/dyac079/6583255?login=true
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' # User could set the next option to number of cores on machine:
#' options(boot.ncpus=2)
#' # Direct and pathway specific attributable fractions estimated
#' # on simulated case control stroke data:
#' # Note that the models here are weighted regressions (based on a column in the
#' # dataframe named 'weights') to rebalance the case control structure to make it
#' # representative over the population, according to the prev argument.
#' # Unweighted regression is fine to use if the data arises from cohort or
#' # cross sectional studies, in which case prev should be set to NULL
#' response_model <- glm(case ~ region * ns(age, df = 5) + sex * ns(age, df = 5) +
#'  education + exercise + ns(diet, df = 3) + smoking + alcohol + stress +
#'   ns(lipids, df = 3) + ns(waist_hip_ratio, df = 3) + high_blood_pressure,
#'   data=stroke_reduced,family='binomial', weights=weights)
#' mediator_models <- list(glm(high_blood_pressure ~ region * ns(age, df = 5) +
#'  sex * ns(age, df = 5) + education   +exercise + ns(diet, df = 3) + smoking +
#'  alcohol + stress,data=stroke_reduced,family='binomial',weights=weights),
#'  lm(lipids ~ region * ns(age, df = 5) + sex * ns(age, df = 5) +education +
#'   exercise + ns(diet, df = 3) + smoking + alcohol + stress, weights=weights,
#'    data=stroke_reduced),lm(waist_hip_ratio ~ region * ns(age, df = 5) +
#'    sex * ns(age, df = 5) + education + exercise + ns(diet, df = 3) +
#'     smoking + alcohol + stress, weights=weights, data=stroke_reduced))
#' # point estimate
#' ps_paf(response_model=response_model, mediator_models=mediator_models ,
#' riskfactor="exercise",refval=0,data=stroke_reduced,prev=0.0035, ci=FALSE)
#' # confidence intervals
#' \donttest{
#' ps_paf(response_model=response_model, mediator_models=mediator_models ,
#' riskfactor="exercise",refval=0,data=stroke_reduced,prev=0.0035, ci=TRUE,
#' boot_rep=100,ci_type="norm")
#'}
ps_paf <- function(response_model, mediator_models,riskfactor,refval,data,prev=NULL,ci=FALSE,boot_rep=100,ci_level=0.95,ci_type=c("norm"), weight_vec=NULL){

  data <- as.data.frame(data)
  N <- nrow(data)
   mediator_names <- c()
  for(i in 1:length(mediator_models)) mediator_names[i] <- as.character(formula(mediator_models[[i]]))[2]
    if(!ci){
      return_vec <- ps_paf_inner(data=data,ind=1:N,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,prev=prev,weight_vec=weight_vec)
      return_vec_names <- c("Direct",mediator_names)
      names(return_vec) <- return_vec_names
      return(return_vec)
    }
  if(ci){
    nc <- options()$boot.ncpus
    cl <- parallel::makeCluster(nc)
    parallel::clusterExport(cl, c("ns"))
    parallel::clusterExport(cl, c("impact_fraction","if_bruzzi","if_direct","predict_df_discrete","pspaf_discrete"))
    res <- boot::boot(data=data,statistic=ps_paf_inner,R=boot_rep,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,prev=prev,weight_vec=weight_vec,cl=cl)
    parallel::stopCluster(cl)
  }
  #options(warn = defaultW)
  return(extract_ci(res=res,model_type='glm',t_vector=c("Direct",mediator_names),ci_level=ci_level,ci_type=ci_type,continuous=TRUE))
}

#' Internal, pathway specific PAF when the mediator is discrete
#'
#' @param data dataframe. A dataframe (with no missing values) containing the data used to fit the mediator and response models.  You can run data_clean to the input dataset if the data has missing values as a pre-processing step
#' @param refval For factor valued risk factors, the reference level of the risk factor.  If the risk factor is numeric, the reference level is assumed to be 0
#' @param riskfactor_col  Integer indicator for the risk factor column in data
#' @param mediator_col Integer indicator for the discrete mediator column in data
#' @param mediator_model A glm or polr model for the mediator, depending on the same confounders and risk factor as specified in the response model.
#' @param response_model A R model object for a binary outcome that involves a risk factor, confounders and mediators of the risk factor outcome relationship.  Note that a weighted model should be used for case control data.  Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param weight_vec A numeric column of weights
#' @return A numeric vector (if ci=FALSE), or data frame (if CI=TRUE) containing estimated PS-PAF for each mediator referred to in mediator_models, together with estimated direct PS-PAF and possibly confidence intervals.
#' @export

pspaf_discrete <- function(data,refval,riskfactor_col,mediator_col,mediator_model,response_model,weight_vec){

  # set up dataframes for prediction (separately for mediator and response)
  inner_bracket <- numeric(nrow(data))
  data_mediator <- data
  if(is.factor(data_mediator[,riskfactor_col])) data_mediator[,riskfactor_col] <- factor(rep(refval,nrow(data)),levels=levels(data_mediator[,riskfactor_col]))
  if(!is.factor(data_mediator[,riskfactor_col])) data_mediator[,riskfactor_col] <- rep(refval,nrow(data))
  if(class(mediator_model)[1]=='polr'){

    mediator_probs <- predict(mediator_model,newdata=data_mediator,type="probs")
  }else{
    # glm

    mediator_probs <- predict(mediator_model,newdata=data_mediator,type="response")
    mediator_probs <- cbind(1-mediator_probs,mediator_probs)

  }
  levels <- sort(unique(data[,mediator_col]))
  if(is.factor(data[,mediator_col])) levels <- levels(data[,mediator_col])
  for(i in 1:length(levels)){
    newd_frame_response <- data
    if(!is.factor(data[,mediator_col])) newd_frame_response[,mediator_col] <- rep(levels[i],nrow(data))  # mediator set to  level i (i could be 1,2,3)
    if(is.factor(data[,mediator_col])) newd_frame_response[,mediator_col] <- factor(rep(levels[i],nrow(data)),levels=levels(data[,mediator_col]))  # mediator set to  level i (i could be 1,2,3)

    predicted_response <- predict(response_model,newdata=newd_frame_response,type="response")
    inner_bracket <- inner_bracket+predicted_response*mediator_probs[,i]

  }

  sum(weight_vec*(predict(response_model,type="response")-inner_bracket))/sum(weight_vec*predict(response_model,type="response"))

}


ps_paf_inner <- function(data, ind, response_model, mediator_models,riskfactor,refval,nsims=1,prev=NULL, weight_vec){


  ###############################
  if(is.null(weight_vec)) weight_vec <- rep(1,nrow(data))
  response <-  as.character(formula(response_model))[2]
  y <- data[,colnames(data)==response]
  N <- nrow(data)
  if(!is.null(prev)){

    data_prev <- mean(as.numeric(y==1))
    weight_vec[y==0] <- (1-prev)/(1-data_prev)
    weight_vec[y==1] <- prev/data_prev

  }
  N <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- class(response_model)[1]

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    mediator_model_type[i] <- class(mediator_models[[i]])[1]
  }
  if(!all(ind==(1:N))){

    weight_vec <- weight_vec[ind]
    data <- data[ind, ]
    y <- data[,colnames(data)==response]
    N <- nrow(data)
    if(!is.null(prev)){

      data_prev <- mean(as.numeric(y==1))
      weight_vec[y==0] <- (1-prev)/(1-data_prev)
      weight_vec[y==1] <- prev/data_prev

    }
    response_model <- update(response_model,data=data)

    for(i in 1:M){
      mediator_models[[i]] <- update(mediator_models[[i]],data=data)

    }

  }
  out_vec <- numeric(M+1)


    new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data=data)
    out_vec[1] <- impact_fraction(model=response_model, data=data, new_data=new_data_direct,calculation_method="D", ci=FALSE, weight_vec=weight_vec)

    for(j in 1:M){

   if(mediator_model_type[j]=='glm' || mediator_model_type[j]=='polr') out_vec[1+j]  <- pspaf_discrete(data=data,refval=refval,riskfactor_col=riskfactor_col,mediator_col=mediator_col[j],mediator_models[[j]],response_model,weight_vec=weight_vec)

   if(mediator_model_type[j]=='lm'){

      mediator_effects <- predict(mediator_models[[j]]) - predict(mediator_models[[j]],new_data_direct)
      new_mediator_data <- data
      new_mediator_data[,mediator_col[j]] <- new_mediator_data[,mediator_col[j]] - mediator_effects
     out_vec[j+1] <- impact_fraction(model=response_model, data=data, new_data=new_mediator_data,calculation_method="D", ci=FALSE, weight_vec=weight_vec)

       }
    }
  return(out_vec)
}






