
ps_paf_sim <- function(response_model, mediator_models,riskfactor,refval,data,prev=NULL,nsims=1,ci=FALSE,boot_rep=100,ci_level=0.95,ci_type=c("norm")){
  N <- nrow(data)
  mediator_names <- c()
  for(i in 1:length(mediator_models)) mediator_names[i] <- as.character(formula(mediator_models[[i]]))[2]
  if(!ci){
    return_vec <- ps_paf_inner_sim(data=data,ind=1:N,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,nsims=nsims,prev=prev)
    return_vec_names <- c("Direct",mediator_names)
    names(return_vec) <- return_vec_names
    return(return_vec)
  }

  if(ci){
    nc <- options()$boot.ncpus
    cl <- parallel::makeCluster(nc)
    parallel::clusterExport(cl, c("ns"))
    res <- boot::boot(data=data,statistic=ps_paf_inner_sim,R=boot_rep,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,nsims=nsims,prev=prev)
    parallel::stopCluster(cl)
  }

  return(extract_ci(res=res,model_type='glm',t_vector=c("Direct",mediator_names),ci_level=ci_level,ci_type=ci_type,continuous=TRUE))
}

#' Estimate Pathway specific population attributable fractions
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
ps_paf <- function(response_model, mediator_models,riskfactor,refval,data,prev=NULL,ci=FALSE,boot_rep=100,ci_level=0.95,ci_type=c("norm")){
  #defaultW <- getOption("warn")
  #options(warn = -1)
  N <- nrow(data)
  mediator_names <- c()
  for(i in 1:length(mediator_models)) mediator_names[i] <- as.character(formula(mediator_models[[i]]))[2]
    if(!ci){
      return_vec <- ps_paf_inner(data=data,ind=1:N,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,prev=prev)
      return_vec_names <- c("Direct",mediator_names)
      names(return_vec) <- return_vec_names
      #options(warn = defaultW)
      return(return_vec)
    }
  if(ci){
    nc <- options()$boot.ncpus
    cl <- parallel::makeCluster(nc)
    parallel::clusterExport(cl, c("ns"))
    res <- boot::boot(data=data,statistic=ps_paf_inner,R=boot_rep,response_model=response_model, mediator_models=mediator_models,riskfactor=riskfactor,refval=refval,prev=prev,cl=cl)
    parallel::stopCluster(cl)
  }
  #options(warn = defaultW)
  return(extract_ci(res=res,model_type='glm',t_vector=c("Direct",mediator_names),ci_level=ci_level,ci_type=ci_type,continuous=TRUE))
}

# function which does one bootstrap rep

ps_paf_inner_sim <- function(data, ind, response_model, mediator_models,riskfactor,refval,nsims=1,prev=NULL){

  #library(splines)
  #library(survival)

  N <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- NULL
  if(grepl("^glm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "glm"
  if(grepl("^lm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "lm"
  if(grepl("^.*polr$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "polr"

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    if(grepl("^glm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "glm"
    if(grepl("^lm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "lm"
    if(grepl("^.*polr$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "polr"
  }

  if(!all(ind==(1:N))){

    data <- data[ind, ]
    if(response_model_type== "glm"){
#browser()
      model_text <- as.character(response_model$call)
      if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"))")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"),weights=",model_text[5],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "lm"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "polr"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }
for(i in 1:M){

if(mediator_model_type[i]== "glm"){
   model_text <- as.character(mediator_models[[i]]$call)
   if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"))")
   if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"),weights=",model_text[5],")")
   model_text <- model_text_u
   thesplit=""
   while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
     model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
     stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
     model_text <- stuff[[1]][1]
     thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
   }
   model_text <- paste0(model_text,thesplit)
   mediator_models[[i]] <- eval(parse(text=model_text_u))
 }

    if(mediator_model_type[i] == "lm"){
      model_text <- as.character(mediator_models[[i]]$call)
      if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      thesplit=""
      while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
        model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
        stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
        model_text <- stuff[[1]][1]
        thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
      }
      model_text <- paste0(model_text,thesplit)
      mediator_models[[i]] <- eval(parse(text=model_text_u))
    }

    if(mediator_model_type[i] == "polr"){
      model_text <- as.character(mediator_models[[i]]$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      thesplit=""
      while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
        model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
        stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
        model_text <- stuff[[1]][1]
        thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
      }
      model_text <- paste0(model_text,thesplit)
      mediator_models[[i]] <- eval(parse(text=model_text_u))
    }
  }



  }
  out_mat <- matrix(0,nrow=M+1,ncol=nsims)

  for(i in 1:nsims){

    new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data=data)
    out_mat[1,i] <- impact_fraction(model=response_model, data=data, new_data=new_data_direct,calculation_method="D", prev=prev,ci=FALSE)

       for(j in 1:M){

      new_data_mediator_j <- data
      new_data_mediator_j[,mediator_col[j]] <- do_sim(colnum=mediator_col[j],current_mat=new_data_direct,model=mediator_models[[j]])
      if(!is.factor(data[,mediator_col[j]]) && is.integer(data[,mediator_col[j]])) new_data_mediator_j[,mediator_col[j]] <- as.integer(new_data_mediator_j[,mediator_col[j]])
      if(is.factor(data[,mediator_col[j]])) new_data_mediator_j[,mediator_col[j]] <- factor(new_data_mediator_j[,mediator_col[j]],levels=levels(data[,mediator_col[j]]))
      out_mat[1+j,i] <- impact_fraction(model=response_model, data=data, new_data=new_data_mediator_j,calculation_method="D", prev=prev,ci=FALSE)

    }

  }
  return(apply(out_mat,1,mean))
}


ps_paf_inner <- function(data, ind, response_model, mediator_models,riskfactor,refval,nsims=1,prev=NULL){


   ########################  load in impact fraction functions:
  impact_fraction <- function(model, data, new_data, calculation_method="B",prev=NULL,ci=FALSE,boot_rep=100,t_vector=NULL, ci_level=0.95, ci_type=c("norm")){

    if(!is.data.frame(data)){
      stop(
        "data must be a dataframe object")
    }

    # remove data not used to fit model
    data <- data[row.names(data) %in% row.names(model.frame(model)),]
    new_data <- new_data[row.names(data) %in% row.names(model.frame(model)),]

    if(!is.data.frame(new_data)){
      stop(
        "new_data must be a dataframe object")
    }

    if(ncol(data)!=ncol(new_data) || nrow(data)!=nrow(new_data)){
      stop(
        "new_data must be the same dimensions as data")
    }

    if(!all(as.character(lapply(data,class))==as.character(lapply(new_data,class))) || !all(colnames(data)==colnames(new_data))){
      stop(
        "Data types and column names in new_data must match data types and column names in data.  To do this, try creating new_data from data")
    }

    if(!calculation_method %in% c("B","D")){
      stop(
        "Calculation of PAF only possible using the (B)ruzzi or (D)irect method.  Please supply either B or D")
    }
    response <- as.character(model$formula)[2]

    model_type <- NULL
    if(grepl("^glm$",as.character(model$call)[1],perl=TRUE)){
      model_type <- "glm"
      if (!as.character(model$family[1])=="binomial" & ! as.character(model$family[2]) %in% c("logit","log")) {
        stop(
          "The family must be binomial and link must be either log or logistic"
        )
      }
    }
    if(grepl("^coxph$",as.character(model$call)[1],perl=TRUE)){
      if("userCall" %in% names(model)){
        model_type <- "clogit"
        vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
        vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^strata\\((.*)\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^Surv\\(rep\\([0-9]*,[0-9]*\\),(.*)\\)$',replacement='\\1',x=vars)
        response <- vars[1]
      }else{
        model_type <- "coxph"
      }
    }
    if (is.null(model_type)) {
      stop(
        "Model must be either a glm, conditional logistic regression or Cox-proportional hazards regression"
      )
    }

    if (model_type=="coxph" && calculation_method=="B") {
      stop(
        "Bruzzi method unavailable with proportional hazards regression due to censoring.  Set method to direct instead"
      )
    }

    if(!is.null(prev) && model_type=="coxph"){

      stop(
        "Prevalence weighted estimation not appropriate for survival data sets"
      )


    }

    if(!is.double(t_vector) && !is.integer(t_vector) && model_type=="coxph"){

      stop(
        "Specify a numeric vector of times at which to calculate PAF"
      )


    }


    if(!is.null(prev) && (prev>=1 || prev<=0)){

      stop(
        "Prevalence should be a proportion (a number between 0 and 1)"
      )
    }

    N <- nrow(data)
    if(calculation_method=="B"){

      if(!ci) return(if_bruzzi(data, ind=1:N, model=model,model_type=model_type,new_data=new_data,response=response))
      if(ci){

        res <- boot::boot(data=data,statistic=if_bruzzi,R=boot_rep, model=model,model_type=model_type,new_data=new_data,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }
    if(calculation_method=="D"){

      if(!ci) return(if_direct(data,ind=1:N,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response))
      if(ci){
        res <- boot::boot(data=data,statistic=if_direct,R=boot_rep,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }

  }



  if_bruzzi <- function(data,ind, model,model_type,  new_data,response){


    N <- nrow(data)

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
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))


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

      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #eta1 <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #eta2 <- the.mat%*%model$coefficients
      #oldRR <- exp(eta1)
      #newRR <- exp(eta2)
      oldRR <- predict(model,type="risk")
      newRR <- predict(model,type="risk",newdata=new_data)
      y <- data[,colnames(data)==response]
    }


    if(model_type == "glm"){

      if(!all(ind==(1:N))){

        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
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

      # predict on linear predictor scale
      oldRR <- exp(predict(model,newdata=data))
      newRR <- exp(predict(model,newdata=new_data))
      y <- data[,colnames(data)==response]
    }

    return(1 - mean(newRR[y==1]/oldRR[y==1]))

  }

  if_direct <- function(data, ind, model,model_type, new_data, prev,t_vector,response){


    N <- nrow(data)
    if(model_type == "coxph"){

      if(!all(ind==(1:N))){
        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("survival::coxph(",model_text[2],",data=data)")
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
      cum_haz <- survival::basehaz(model, centered=FALSE)
      t_indices <- integer(length(t_vector))
      for(i in 1:length(t_vector)){
        t_indices[i] <- which.min(sapply(cum_haz[,2],function(x){(x-t_vector[i])^2}))
      }
      cum_haz <- cum_haz[t_indices,]
      oldhr <- predict(model,type="risk")
      newhr <- predict(model,newdata=new_data,type="risk")

      mean_probs_old <- 1 - apply(exp(-outer(cum_haz[,1],oldhr)),1,mean)

      mean_probs_new <- 1 - apply(exp(-outer(cum_haz[,1],newhr)),1,mean)

      PAF_vec <- (mean_probs_old - mean_probs_new)/mean_probs_old
      names(PAF_vec) <- paste0("t=",round(cum_haz[,2],2))
      return(PAF_vec)

    }



    add_term <- 0

    if(model_type=="clogit"){

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
        new_data <- new_data[resamples,]
        data <- data[resamples,]
        # avoid duplication of strata names
        data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))

        #refit model

        model_text <- paste0("survival::clogit(",model_text,",data=data)")
        thesplit=""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
        model <- eval(parse(text=model_text))

      }
      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #lp_old <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #lp_new <- the.mat%*%model$coefficients
      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model, newdata=new_data)
      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(!is.null(prev)){

        temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
        add_term <- uniroot(temp_fun, interval=c(-100,100))$root

      }
      probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
      probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))


    }else{
      # model is a glm

      if(!all(ind==(1:N))){
        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
        model_text <- model_text_u
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

      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model,newdata=new_data)

      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(as.character(model$family[2])=="logit"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
        probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))

      }
      if(as.character(model$family[2])=="log"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)
        probs_new <- exp(lp_new+add_term)

      }
    }
    return((sum(weights*probs_old)-sum(weights*probs_new))/sum(weights*probs_old))

  }

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


  pspaf_discrete <- function(data,refval,riskfactor_col,mediator_col,mediator_model,response_model,weights){

    # set up dataframes for prediction (separately for mediator and response)
    inner_bracket <- numeric(nrow(data))
    data_mediator <- data
    if(is.factor(data_mediator[,riskfactor_col])) data_mediator[,riskfactor_col] <- factor(rep(refval,nrow(data)),levels=levels(data_mediator[,riskfactor_col]))
    if(!is.factor(data_mediator[,riskfactor_col])) data_mediator[,riskfactor_col] <- rep(refval,nrow(data))
    if(names(mediator_model)[2]=='zeta'){

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

    sum(weights*(predict(response_model,type="response")-inner_bracket))/sum(weights*predict(response_model,type="response"))

  }


  ###############################

  N <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- NULL
  if(grepl("^glm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "glm"
  if(grepl("^lm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "lm"
  if(grepl("^.*polr$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "polr"

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    if(grepl("^glm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "glm"
    if(grepl("^lm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "lm"
    if(grepl("^.*polr$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "polr"
  }

  if(!all(ind==(1:N))){

    data <- data[ind, ]
    if(response_model_type== "glm"){
      #browser()
      model_text <- as.character(response_model$call)
      if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"))")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"),weights=",model_text[5],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "lm"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "polr"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }
    for(i in 1:M){

      if(mediator_model_type[i]== "glm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"),weights=",model_text[5],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "lm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "polr"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }
    }



  }
  out_vec <- numeric(M+1)


    new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data=data)
    out_vec[1] <- impact_fraction(model=response_model, data=data, new_data=new_data_direct,calculation_method="D", prev=prev,ci=FALSE)

    for(j in 1:M){

   if(mediator_model_type[j]=='glm' || mediator_model_type[j]=='polr') out_vec[1+j]  <- pspaf_discrete(data=data,refval=refval,riskfactor_col=riskfactor_col,mediator_col=mediator_col[j],mediator_models[[j]],response_model,weights=mediator_models[[j]]$weights)

   if(mediator_model_type[j]=='lm'){

      mediator_effects <- predict(mediator_models[[j]]) - predict(mediator_models[[j]],new_data_direct)
      new_mediator_data <- data
      new_mediator_data[,mediator_col[j]] <- new_mediator_data[,mediator_col[j]] - mediator_effects
     out_vec[j+1] <- impact_fraction(model=response_model, data=data, new_data=new_mediator_data,calculation_method="D", prev=prev,ci=FALSE)

       }
    }
  return(out_vec)
}






