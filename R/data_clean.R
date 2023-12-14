#' Clean a dataset to make model fitting more efficient
#'
#' Strip out unneeded variables from original data (based on fitted model, or alternatively based on specifying a list of variables), and remove rows with NA values.  The function works for logistic, survival and conditional logistic regressions.  The function also creates a column of weights, which will be just a vector of 1s if prevalence is unspecified.
#'
#' @param model A glm (with logistic or log link, with binomial family), clogit or coxph model.
#' @param data A data frame that was used to fit the model
#' @param vars Default NULL.  Variables required in output data set.  If set to NULL and model is specified, the variables kept are the response and covariates assumed in model.  If set to NULL and model is unspecified, the original dataset is returned.
#' @param response Default "case".  response variable in dataset.  Used when recalculating weights (if the argument prev is set)  If set to NULL, the response is inferred from the model
#' @param prev Default NULL.  Prevalence of disease (or yearly incidence of disease in healthy controls).  Only relevant to set in case control studies and if path specific PAF or sequential joint PAF calculations are required.  The purpose of this is to create a vector of weights in output dataset, that reweights the cases and controls to reflect the general population.  This vector of weights can be used to fit weighted regression models.
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
data <- as.data.frame(data)
if(is.null(vars)){
  if(is.null(model)) return(data)
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


