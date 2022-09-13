#' General calculations of impact fractions
#'
#' @param model Either a clogit, glm or survival regression
#' @param data A dataframe containing variables used for fitting the model
#' @param new_data A dataframe (of the same variables and size as data) representing an alternative distribution of risk factors
#' @param calculation_method A character either 'B' (Bruzzi) or 'D' (Direct method).  For case control data, the method described in Bruzzi 1985 is recommended.  Bruzzi's method estimates PAF from relative risks and prevalence of exposure to the risk factor.  The Direct method estimates PAF by summing estimated probabilities of disease in the absense of exposure on the individual level
#' @param prev estimated prevalence of disease.  This only needs to be specified if the data source is from a case control study, and the direct method is used
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with point estimate (default FALSE)
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE)
#' @param t_vector Numeric.  A vector of times at which to calculate PAF (only specified if model is coxph)
#' @param ci_level Numeric.  Default 0.95. A number between 0 and 1 specifying the confidence level
#' @param ci_type Character.  Defalt norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @return An estimated impact if ci=FALSE, or for survival data a vector of estimated impact corresponding to event times in the data.  If ci=TRUE, a vector with elements corresponding to the raw estimated impact fraction, estiamted bias, bias corrected estimate and lower and upper elements of any confidence procedures requested.  If ci=TRUE, and a coxph model is fit, a matrix will be returned, with rows corresponding to differing times at which the impact fraction might be calculated.
#' @export
#'
#' @examples
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

  library(splines)
  library(survival)

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
      model_text <- paste0("clogit(",model_text,",data=data)")
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

  library(splines)
  library(survival)

  N <- nrow(data)
  if(model_type == "coxph"){

    if(!all(ind==(1:N))){

      data <- data[ind, ]
      new_data <- new_data[ind, ]
      model_text <- as.character(model$call)
      model_text <- paste0("coxph(",model_text[2],",data=data)")
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
      model_text <- paste0("clogit(",model_text,",data=data)")
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
      model <- eval(parse(text=model_text_u))

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

extract_ci <- function(res,model_type,t_vector,ci_level,ci_type,continuous=FALSE){
if(continuous){

  d <- data.frame(matrix(ncol=3 + 2*length(ci_type),nrow=length(res$t0)))
  colnames(d) <- c("raw_estimate", "bias","bias_corrected",rep("",2*length(ci_type)))
  for(i in 1:length(ci_type)) colnames(d)[(2+2*i):(3+2*i)] <- c(paste0(ci_type[i],"_lower"),paste0(ci_type[i],"_upper"))



  d[,1] <- res$t0
  d[,2] <- apply(res$t,2,mean,na.rm=TRUE)-res$t0
  d[,3] <- 2*res$t0-apply(res$t,2,mean,na.rm=TRUE)
  for(j in 1:length(res$t0)){
    for(i in 1:length(ci_type)){

      if(ci_type[i]=="norm") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="norm", index=j)$normal[2:3]
      if(ci_type[i]=="basic") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="basic", index=j)$basic[4:5]
      if(ci_type[i]=="perc") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="perc", index=j)$perc[4:5]
      if(ci_type[i]=="bca") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="bca", index=j)$bca[4:5]
    }
  }
  rownames(d) <- t_vector
  for(i in 1:ncol(d)) d[,i] <- signif(d[,i],3)
  return(d)

}
  if(model_type!="coxph"){

    v <- numeric(3 + 2*length(ci_type))
    names(v) <- c("raw_estimate", "estimated_bias","bias_corrected_estimate",rep("",2*length(ci_type)))
    for(i in 1:length(ci_type)) names(v)[(2+2*i):(3+2*i)] <- c(paste0(ci_type[i],"_lower"),paste0(ci_type[i],"_upper"))
    v[1] <- res$t0
    v[2] <- mean(res$t)-res$t0
    v[3] <- 2*res$t0-mean(res$t)

    for(i in 1:length(ci_type)){
      if(ci_type[i]=="norm") v[(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="norm")$normal[2:3]
      if(ci_type[i]=="basic") v[(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="basic")$basic[4:5]
      if(ci_type[i]=="perc") v[(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="perc")$perc[4:5]
          if(ci_type[i]=="bca") v[(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="bca")$bca[4:5]
    }
    v <- signif(v,digits=3)
    return(v)

  }
  if(model_type=="coxph"){

    d <- data.frame(matrix(ncol=3 + 2*length(ci_type),nrow=length(t_vector)))
    colnames(d) <- c("raw_estimate", "estimated_bias","bias_corrected_estimate",rep("",2*length(ci_type)))
    for(i in 1:length(ci_type)) colnames(d)[(2+2*i):(3+2*i)] <- c(paste0(ci_type[i],"_lower"),paste0(ci_type[i],"_upper"))



    d[,1] <- res$t0
    d[,2] <- apply(res$t,2,mean)-res$t0
    d[,3] <- 2*res$t0-apply(res$t,2,mean)
    for(j in 1:length(t_vector)){
    for(i in 1:length(ci_type)){

      if(ci_type[i]=="norm") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="norm", index=j)$normal[2:3]
      if(ci_type[i]=="basic") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="basic", index=j)$basic[4:5]
      if(ci_type[i]=="perc") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="perc", index=j)$perc[4:5]
      if(ci_type[i]=="bca") d[j,(2+2*i):(3+2*i)] <- boot::boot.ci(res, conf=ci_level,type="bca", index=j)$bca[4:5]
    }
    }
    rownames(d) <- t_vector
    for(i in 1:ncol(d)) d[,i] <- signif(d[,i],3)
    return(d)
  }
}
