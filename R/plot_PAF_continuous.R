#' Plot impact fractions corresponding to risk-quantiles over several risk factors
#'
#' @param x A PAF_q object.  This is a dataframe that is created by running the function PAF_calc_continuous.
#' @param ... Other global arguments inherited by that might be passed to the ggplot routine
#' @return A ggplot2 plotting object for PAF_q over the differing risk factors in x
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
#'   high_blood_pressure, family = "binomial", data = stroke_reduced)
#' out <- PAF_calc_continuous(model_continuous,riskfactor_vec=
#' c("diet","lipids","waist_hip_ratio"),q_vec=c(0.1,0.9),
#' ci=FALSE,calculation_method="B",data=stroke_reduced)
#' plot(out)
#' \donttest{
#' # example with more quantile points and including confidence intervals
#' # (more useful - but a bit slower to run)
#' out <- PAF_calc_continuous(model_continuous,riskfactor_vec=
#' c("diet","lipids","waist_hip_ratio"),q_vec=c(0.01, 0.1,0.3,0.5,0.7,0.9),
#' ci=TRUE,calculation_method="B",data=stroke_reduced)
#' plot(out)
#' }
plot.PAF_q <- function(x, ...){
   #options(warn = -1)

  if(!"riskfactor" %in% names(x)){

    stop("data_frame should have a column named 'riskfactor' with the names of the risk factors or exposures being compared")

  }

  if(!"q_val" %in% names(x)){

    stop("data_frame should have a column named 'q_val' specifying the intervention on the risk factor that is considered")

  }
  data_frame <- data.frame(riskfactor=x$riskfactor,q = x$q_val, paf_q=signif(x$paf_q,3))
  riskfactor <- x$riskfactor
  q_val <- x$q_val

  pd <- ggplot2::position_dodge(0.01)
  J <- ncol(data_frame)
  if(J > 3){
  colnames(data_frame)[(J-2):J] <- c("PAF","lower_CI", "upper_CI")
p <- ggplot2::ggplot(data_frame, ggplot2::aes(x=-q_val, y=PAF, colour=riskfactor)) + ggplot2::theme_grey(base_size = 20) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=lower_CI, ymax=upper_CI),width=0.1, position=pd,size=1) + ggplot2::geom_line(ggplot2::aes(group=riskfactor,width=1),position=pd,size=1) +
  ggplot2::labs(col = "Risk Factor")+ggplot2::ylab(expression(paste(hat(PAF[q]))))+ggplot2::xlab("q") + ggplot2::scale_x_continuous(breaks=c(-1,-.8,-0.6,-0.4,-0.2,0),labels=c(1,0.8,0.6,0.4,0.2,0))

  }
  if(J==3){
    colnames(data_frame)[J] <- "PAF"
   p <-  ggplot2::ggplot(data_frame, ggplot2::aes(x=-q_val, y=PAF, colour=riskfactor)) + ggplot2::theme_grey(base_size = 20) + ggplot2::geom_line(ggplot2::aes(group=riskfactor,width=1),position=pd,size=1) +
      ggplot2::labs(col = "Risk Factor")+ggplot2::ylab(expression(paste(hat(PAF[q]))))+ggplot2::xlab("q") + ggplot2::scale_x_continuous(breaks=c(-1,-.8,-0.6,-0.4,-0.2,0),labels=c(1,0.8,0.6,0.4,0.2,0))

  }
   p
}



plot_continuous_quick <- function(model,riskfactor,data,S = 10,ref_val=NA, ci_level=0.95,min_risk_q=.1,n_x=1000,plot_region=TRUE, plot_density=TRUE,theylab="OR",qlist){
   data_orig <- data
  if(n_x<=nrow(data)){
    data <- data[sample(1:nrow(data),n_x,replace=FALSE),]
  }
  xmat <- model.matrix(model, data=data) # design matrix
  col_indexes <- grep(paste('^.*',riskfactor,'.*$',sep=''),
                      colnames(xmat),perl=TRUE)
  # calculate submatrix necessary for plot
  submat <- xmat[,col_indexes]
  # values of x in data
  riskfactor_col <- grep(paste("^",riskfactor,"$",sep=''),colnames(data),perl=TRUE)
  xvals <- data[,riskfactor_col]

  model_frame <- model.matrix(model,data=data)

  model_frame <- as.matrix(model_frame)

  M <- median(xvals)
  index <- which.min((xvals-M)^2)
  if(!is.na(ref_val)){

    index <- which.min((xvals-ref_val)^2)
  }

  median_mat <- matrix(rep(model_frame[index,],nrow(data)),nrow=nrow(data),byrow=TRUE)
  model_frame <- model_frame-median_mat
  pred <- model_frame[,col_indexes,drop=FALSE]%*%model$coefficients[col_indexes,drop=FALSE]
  var <- diag(model_frame[,col_indexes,drop=FALSE]%*%vcov(model)[col_indexes,col_indexes,drop=FALSE]%*%t(model_frame[,col_indexes,drop=FALSE]))
  z_val <- qnorm(1-(1-ci_level)/2)
  lower <- exp(pred - z_val*sqrt(var))
  upper <- exp(pred + z_val*sqrt(var))
  middle <- exp(pred)
  d <- data.frame(x=xvals,est=middle,lower=lower,upper=upper)
  d <- d[order(d$x,decreasing=FALSE),]
  gg <- ggplot2::ggplot(d, ggplot2::aes(x=x, y=middle)) + ggplot2::geom_line(ggplot2::aes(x=x, y=est),lwd=2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),alpha=0.5,fill='red')+ ggplot2::theme_grey(base_size = 20) +ggplot2::xlab(riskfactor)+ggplot2::ylab(theylab)


  if(plot_region){


    risk_quantiles <- risk_quantiles(riskfactor = riskfactor,data=data_orig,model=model,S=S,q=qlist)

    d$min=0

    d$max=0
    for(j in 1:nrow(d)){
      d$max[j] <- d$est[j]
      index <- which.min((d$x[j]-risk_quantiles)^2)
      if(as.numeric(names(risk_quantiles)[index])>min_risk_q) d$max[j] <- 0
    }
    #polygon(themat, col=adjustcolor( "blue", alpha.f = .2))
    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=min, ymax=max),alpha=0.5,fill='blue',data=d)
  }
  if(plot_density){
    d1 <- density(xvals)
    scaling <- max(d$est)*0.7/max(pmax(d1$y))
    xdensity=d1$x
    ydensity=scaling*d1$y
    newd <- data.frame(xdensity=xdensity,ydensity=ydensity)
    newd <- newd[xdensity>=quantile(xvals,min(qlist))&xdensity<=quantile(xvals,max(qlist)),]
    gg <- gg + ggplot2::geom_line(ggplot2::aes(x=xdensity,y=ydensity),lwd=2,col='black',data=newd)
  }
   gg

}


#' Plot hazard ratios, odds ratios or risk ratios comparing differing values of a continuous exposure to a reference level
#' @param model A fitted model (either glm, clogit or coxph)
#' @param riskfactor The string name of a continuous exposure or risk factor represented in the data and model
#' @param data Data frame used to fit the model
#' @param ref_val The reference value used in plotting.  If left at NA, the median value of the risk factor is used
#' @param ci_level Numeric.  A number between 0 and 1 specifying the confidence level
#' @param min_risk_q Default .1.  A number between 0 and 1 representing the desired risk quantile for the continuous risk factor
#' @param plot_region Default TRUE.  Logical specifying whether the targeted region corresponding to an intervention setting the continuous risk factor at a quantile min_risk_q or lower is to be plotted
#' @param plot_density  Default TRUE.  Logical specifying whether density of distribution of risk factor is to be added to the plot
#' @param n_x  Default 10000.  How many values of riskfactor will be used to plot spline (when interact=FALSE)
#' @param theylab  Default "OR".  Y-axis label of the plot
#' @param qlist Vector of quantile values for q, corresponding to the plotted values of PAF_q for each risk factor/exposure
#' @param interact  Default "FALSE".  Set to TRUE spline models enter as interactions.
#' @param S  Default 10.  The integer number of random samples used to calculate average differences in linear predictors. Only relevant to set when interact=TRUE
#' @return A ggplot2 plotting object
#' @export
#' @references Ferguson, J., Maturo, F., Yusuf, S. and O’Donnell, M., 2020. Population attributable fractions for continuously distributed exposures. Epidemiologic Methods, 9(1)
#'
#' @examples
#' library(survival)
#' library(splines)
#' model_continuous <- glm(formula = case ~ region * ns(age, df = 5) +
#'  sex * ns(age, df = 5) + education +exercise + ns(diet, df = 3) +
#'  alcohol + stress + ns(lipids,df = 3) + ns(waist_hip_ratio, df = 3) +
#'   high_blood_pressure, family = "binomial", data = stroke_reduced)
#' plot_continuous(model_continuous,riskfactor="diet",data=stroke_reduced)
plot_continuous <- function(model,riskfactor,data,S = 10,ref_val=NA, ci_level=0.95,min_risk_q=.1,plot_region=TRUE, plot_density=TRUE,n_x=10000,theylab="OR",  qlist=seq(from=0.001,to=0.999,by=0.001), interact=FALSE){

  if(!interact){
        return(plot_continuous_quick(model=model,riskfactor=riskfactor,data=data,ci_level=ci_level,min_risk_q=min_risk_q,n_x=n_x,plot_region=plot_region,plot_density=plot_density,theylab=theylab,qlist=qlist))
  }
  #options(warn = defaultW)
  #options(warn = -1)
  xmat <- model.matrix(model, data=data) # design matrix
  col_indexes <- grep(paste('^.*',riskfactor,'.*$',sep=''),
                      colnames(xmat),perl=TRUE)
  # calculate submatrix necessary for plot
  submat <- xmat[,col_indexes]
  # values of x in data
  riskfactor_col <- grep(paste("^",riskfactor,"$",sep=''),colnames(data),perl=TRUE)
  xvals <- data[,riskfactor_col]


  # plot will be made at .001,.002,...,.999 quantiles of risk factor
  plot_x_v <- quantile(xvals,qlist,na.rm=T)

  # if sampling make dataframe same size as original by default
    newS <- S
    if(S<2) newS <- 2  # model.frame won't work on a dataset of size 1
    data_s <- data[sample(1:nrow(data),newS),,drop=FALSE]
    ind <- rep(1,S)
    for(j in 1:ncol(data_s)){

      if(is.factor(data[,j])) data_s[,j] <- factor(data_s[,j],levels=levels(data[,j]))
      if(is.numeric(data[,j])) data_s[,j] <- as.numeric(data_s[,j])
    }
    newd <- data_s
    newd[,riskfactor_col] <- plot_x_v[1]
    a <- try(model_frame <-model.matrix(model,data=newd),silent=TRUE)
 if(!substring(a[1],1,5)=="Error"){
   model_frame <- model.matrix(model,data=newd)
    for(i in 2:length(qlist)){
      newd <- data_s
      newd[,riskfactor_col] <- plot_x_v[i]
      model_frame <- rbind(model_frame,model.matrix(model,data=newd))
      ind <- c(ind, rep(i,S))
    }
}else
  {
      N <- nrow(data)
      model_frame <- model.matrix(model,data=rbind(newd,data[1:N,]))
      model_frame <- model_frame[(1:(nrow(model_frame)-N)),]
      for(i in 2:length(qlist)){
        newd <- data_s
        newd[,riskfactor_col] <- plot_x_v[i]
        model_frame <- rbind(model_frame,model.matrix(model,data=rbind(newd,data[1:N,])))
        model_frame <- model_frame[(1:(nrow(model_frame)-N)),]
        ind <- c(ind, rep(i,S))

      }
    }

  # collapse model_frame now over ind
  model_frame <- as.data.frame(model_frame)
  model_frame$ind <- ind
  model_frame <- dplyr::group_by(model_frame,ind)
  stuff <- dplyr::funs(mean)
  model_frame <- dplyr::summarize_all(model_frame,stuff)
  model_frame <- model_frame[,colnames(model_frame)!="ind"]
  model_frame <- as.matrix(model_frame)
  # new predictions (pred- pred at median)
  i1 <- length(unique(ind))/2

  if(i1%%1 == 0) i2 <- i1+1
  if(i1%%1 == 0.5){
    i1 <- i1+0.5
    i2 <- i1-11
  }
  median <- (0.5*model_frame[i1,,drop=FALSE] + 0.5*model_frame[i1,,drop=FALSE])
  if(!is.na(ref_val)){

    if(sum(plot_x_v<ref_val)==0 | sum(plot_x_v<ref_val)>=length(plot_x_v)){

      stop(paste0("specify reference riskfactor between ", signif(plot_x_v[1],3)," and ", signif(plot_x_v[length(plot_x_v)],3) ))

    }

    l <- length(plot_x_v[plot_x_v<ref_val])
    w <- (ref_val-plot_x_v[l])/(plot_x_v[l+1]-plot_x_v[l])
    median <- w*model_frame[l+1,,drop=FALSE]+(1-w)*model_frame[l,,drop=FALSE]

  }
  median_mat <- matrix(rep(median,length(qlist)),nrow=length(qlist),byrow=TRUE)
  model_frame <- model_frame-median_mat
  pred <- model_frame[,col_indexes,drop=FALSE]%*%model$coefficients[col_indexes,drop=FALSE]
  var <- diag(model_frame[,col_indexes,drop=FALSE]%*%vcov(model)[col_indexes,col_indexes,drop=FALSE]%*%t(model_frame[,col_indexes,drop=FALSE]))
  z_val <- qnorm(1-(1-ci_level)/2)
  lower <- exp(pred - z_val*sqrt(var))
  upper <- exp(pred + z_val*sqrt(var))
  middle <- exp(pred)
  d <- data.frame(x=plot_x_v,est=middle,lower=lower,upper=upper)
  gg <- ggplot2::ggplot(d, ggplot2::aes(x=plot_x_v, y=middle)) + ggplot2::geom_line(ggplot2::aes(x=plot_x_v, y=middle),lwd=2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),alpha=0.5,fill='red')+ ggplot2::theme_grey(base_size = 20) +ggplot2::xlab(riskfactor)+ggplot2::ylab(theylab)


  if(plot_region){


  risk_quantiles <- risk_quantiles(riskfactor = riskfactor,data=data,model=model,S=S,q=qlist)

  d$min=0

  d$max=0
  for(j in 1:nrow(d)){
    d$max[j] <- d$est[j]
    index <- which.min((d$x[j]-risk_quantiles)^2)
    if(as.numeric(names(risk_quantiles)[index])>min_risk_q) d$max[j] <- 0
  }
  #polygon(themat, col=adjustcolor( "blue", alpha.f = .2))
  gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=min, ymax=max),alpha=0.5,fill='blue',data=d)
  }
  if(plot_density){
    d1 <- density(xvals)
  scaling <- max(d$est)*0.7/max(pmax(d1$y))
  xdensity=d1$x
  ydensity=scaling*d1$y
  newd <- data.frame(xdensity=xdensity,ydensity=ydensity)
  newd <- newd[xdensity>=quantile(xvals,min(qlist))&xdensity<=quantile(xvals,max(qlist)),]
  gg <- gg + ggplot2::geom_line(ggplot2::aes(x=xdensity,y=ydensity),lwd=2,col='black',data=newd)
  }
  #options(warn = defaultW)
  gg

}

