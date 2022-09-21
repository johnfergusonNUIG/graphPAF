#' Clean a dataset to make modeling fitting more efficient
#'
#' Strip out unneeded variables from original data (based on fitted model, or alternatively based on specifying a list of variables), and remove rows with NA values.  The function works for logistic, survival and conditional logistic regressions.  The function also creates a column of weights, which will be just a vector of 1s if prevalence is unspecified.
#'
#' @param model A glm (with logistic or log link, with binomial family), clogit or coxph model.
#' @param data A data frame that was used to fit the model
#' @param vars Default NULL.  Variables required in output data set.  If set to NULL and model is specified, the variables kept are the response and covariates assumed in model
#' @param response Default "case".  response variable in dataset.  Used when recalculating weights (if the argument prev is set)  If set to NULL, the response is inferred from the model
#' @param prev Default NULL.  Prevalence of disease (or yearly incidence of disease in healthy controls).  Only relevant to set in case control studies and if path specific PAF or sequential joint PAF calculations are required.  The purpose of this is to create a vector of weights that reweights the cases and controls to reflect the general population
#' @return A cleaned data frame
#' @export
#' @examples
#' # example of using dataclean to strip out NAs, redundant columns and recalculate weights
#' library(survival)
#' library(splines)
#' stroke_reduced_2 <- stroke_reduced
#' stroke_reduced_2$case[sample(1:length(stroke_reduced_2$case),50)] <- NA
#' stroke_reduced_2$random <- rnorm(length(stroke_reduced_2$case))
#' stroke_reduced_3 <- data_clean(stroke_reduced_2,vars=colnames(stroke_reduced),prev=0.01)
#' dim(stroke_reduced_2)
#' dim(stroke_reduced_3)
#' mymod <- clogit(case ~ high_blood_pressure + strata(strata),data=stroke_reduced_2)
#' stroke_reduced_3 <- data_clean(stroke_reduced_2,model=mymod,prev=0.01)
#' dim(stroke_reduced_2)
#' dim(stroke_reduced_3)
data_clean <- function(data,model=NULL,vars=NULL,response="case", prev=NULL){
if(is.null(vars)){
  model_type <- NULL
  vars <- c()
  if(grepl("^glm$",as.character(model$call)[1],perl=TRUE)){
    model_type <- "glm"
    if (!as.character(model$family[1])=="binomial" & ! as.character(model$family[2]) %in% c("logit","log")) {
      stop(
        "The family must be binomial and link must be either log or logistic"
      )
    }
      vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
      vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
      if(length(as.character(model$call))==5){
          colnames(data)[colnames(data)==as.character(model$call)[5]] <- "weights"
          vars <- c(vars, "weights")
      }
          vars <- unique(vars)
    }

  if(grepl("^coxph$",as.character(model$call)[1],perl=TRUE)){
    if("userCall" %in% names(model)){
      model_type <- "clogit"

      vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
      vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^strata\\((.*)\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^Surv\\(rep\\([0-9]*,[0-9]*\\),(.*)\\)$',replacement='\\1',x=vars)
      vars <- unique(vars)
    }else{
      model_type <- "coxph"
      vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
      vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^strata\\((.*)\\)$',replacement='\\1',x=vars)
      vars <- gsub(pattern='^Surv\\((.*),(.*)\\)$',replacement='\\1,\\2',x=vars)
      vars <- c(unlist(strsplit(vars[1],split="[,]")),vars[2:length(vars)])
      vars <- unique(vars)
    }
  }
}
  if("weights"%in%colnames(data) && !("weights" %in% vars)) vars <- c(vars, "weights")
  data <- data[,colnames(data)%in%vars]
  tokeep <- apply(data,1,function(x){sum(is.na(x))==0})
  data <- data[tokeep,]
   if(!is.null(model)) y <- model$y
   if(!is.null(response)) y <- data[,colnames(data)==response]
    N <- nrow(data)

    if(!is.null(prev)){

      weights <- numeric(nrow(data))
      data_prev <- mean(y)
      weights[y==0] <- (1-prev)/(1-data_prev)
      weights[y==1] <- prev/data_prev
    data$weights <- weights

    }
    if(!(c("weights") %in% colnames(data))) data$weights <- rep(1, N)
  ##  make 2 level factor variables into 0/1 variables (0 being reference level), 2 level numeric variables are converted to 0/1 (with the lower value being reference), 2 level character variables are converted to 0/1 with the value appearing in the data first being the refernece.  Reference is first alphabetical level for character variables.
S <- ncol(data)-1
  for(i in 1:S){
    if(is.factor(data[,i]) && length(levels(data[,i]))==2) data[,i] <- factor(as.numeric(data[,i]==levels(data[,i])[2]),levels=c(0,1))
       if(!colnames(data)[i]=="weights" && !is.factor(data[,i]) && length(unique(data[,i]))==2 && is.numeric(data[,i])) data[,i] <- factor(as.numeric(data[,i]==max(data[,i])), levels=c(0,1))
    if(is.character(data[,i])) data[,i] <- factor(data[,i],levels=sort(unique(data[,i])))

  }
  data
}


refit <- function(model,data,with_weights=FALSE){
  model_type <- NULL
  if(grepl("^glm$",as.character(model$call)[1],perl=TRUE)) model_type <- "glm"
  if(grepl("^lm$",as.character(model$call)[1],perl=TRUE)) model_type <- "lm"
  if(grepl("^.*polr$",as.character(model$call)[1],perl=TRUE)) model_type <- "polr"
  if(grepl("^coxph$",as.character(model$call)[1],perl=TRUE)){
    if("userCall" %in% names(model)){
      model_type <- "clogit"
    }else{
      model_type <- "coxph"
    }
  }
  if(model_type=="clogit"){
    model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
    model_text <- paste0(model_text[2],model_text[1],model_text[3])
    model_text <- paste0("clogit(",model_text,",data=data)")
    model <- eval(parse(text=model_text))
  }
  if(model_type=="coxph"){

    model_text <- as.character(model$call)
    model_text <- paste0("coxph(",model_text[2],",data=data)")
    model <- eval(parse(text=model_text))
  }

  if(model_type== "glm"){
    #browser()
    model_text <- as.character(model$call)
    if(with_weights==FALSE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
    if(with_weights==TRUE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=weights)")
    if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
    model <- eval(parse(text=model_text_u))
  }

  if(model_type == "lm"){
    model_text <- as.character(model$call)
    if(with_weights==FALSE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
    if(with_weights==TRUE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data,weights=weights)")
        if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
    model <- eval(parse(text=model_text_u))
  }

  if(model_type == "polr"){
    model_text <- as.character(model$call)
    if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
    if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
    model <- eval(parse(text=model_text_u))
  }
model
}
