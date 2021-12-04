
#' Create a summary data frame for risk factors
#'
#' Create a summary data frame for risk factors, prevalence and risk ratios.  This will be used in fan plots and nomograms
#'
#' @param rf_names character
#' @param rf_prev double
#' @param risk double
#' @param log logical
#'
#' @return A rf.data.frame object
#' @export
#'
#' @examples
#'
#' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA','Diet','WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)


rf_summary <- function(rf_names, rf_prev, risk, log=TRUE){
  stopifnot(length(rf_names)==length(rf_prev) & length(rf_prev)==length(risk))
  stopifnot(is.character(rf_names))
  stopifnot(all(rf_prev > 0 & rf_prev <1))
  stopifnot(is.double(risk))
  stopifnot((risk>0) | log==TRUE)
  log_riskratio <- risk
  if(!log) log_riskratio <- log(risk)
   approx_PAF <- rf_prev*log_riskratio
  rf_summary <- data.frame(rf_prev=rf_prev,log_riskratio=log_riskratio, approx_PAF=approx_PAF,row.names=rf_names)
  rf_summary <- rf_summary[order(rf_summary$approx_PAF,decreasing=TRUE),]
  N <- length(row.names(rf_summary))
  rf_summary$rf_names <- factor(1:N,labels=paste(1:N,": ", row.names(rf_summary),sep=''))
  structure(as.list(rf_summary),class="rf.data.frame",row.names=row.names(rf_summary))
}


#' Create a fan_plot of a rf.data.frame object
#'
#' Create a fan plot displaying approximate PAF, risk factor prevalence and risk ratios
#'
#' @param rf_data_frame A rf.data.frame object
#'
#' @return
#' @export
#'
#' @examples
#'
#' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA','Diet','WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)
#' fan_plot(rfs)

fan_plot <- function(rf_data_frame){

    if(!class(rf_data_frame)=='rf.data.frame'){

    stop("Create a valid rf.data.frame object before running function")

  }
  rf_data_frame <- structure(as.list(rf_data_frame),class="data.frame", row.names=attr(rf_data_frame,"row.names"))

  rf_data_frame$inv_prev <- 1/rf_data_frame$rf_prev
  p <- ggplot2::ggplot(rf_data_frame, ggplot2::aes(inv_prev,log_riskratio)) + ggplot2::geom_point(size=4,ggplot2::aes(color=rf_data_frame$rf_names))  # X and Y axis limits

  p <- p + ggplot2::labs(col = "Risk factors (ranked)")
  for(i in 1:dim(rf_data_frame)[1]){

    temprf_data_frame <- data.frame(x=c(1,rf_data_frame$inv_prev[i]),y=c(rf_data_frame$approx_PAF[i],rf_data_frame$log_riskratio[i]))
    p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y), linetype="dashed", size=1.5,col="black")
  }

  p <- p + ggplot2::theme_minimal()
## ad scale here
  p <- p + ggrepel::geom_label_repel(ggplot2::aes(label=rf_data_frame$rf_names,x=inv_prev+.5,y=log_riskratio), size=4, data=rf_data_frame)

  p <- p + ggplot2::scale_x_continuous(breaks=c(1,1/0.5,1/.3,1/.2,1/.1,1/.05), minor_breaks=NULL,labels=c("100%","50%","30%","20%","10%","5%"),name="Prevalence in Controls")


  p <- p + ggplot2::scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1),labels=c(paste0(seq(from=0,to=100,by=10),"%"),""),minor_breaks=NULL,sec.axis=ggplot2::sec_axis(~.,breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1),labels=floor(exp(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1))*10)/10,name="Odds Ratio"),name="approximate PAF")

  p <- p+ ggplot2::theme(axis.text.x = ggplot2::element_text(colour="grey20",size=30,angle=90,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = ggplot2::element_text(colour="blue",size=30,angle=0,hjust=1,vjust=0,face="plain"),        axis.text.y.right = ggplot2::element_text(colour="grey20",size=30,angle=0,hjust=1,vjust=0,face="plain"),
          axis.title.x = ggplot2::element_text(colour="grey20",size=30,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = ggplot2::element_text(colour="blue",size=30,angle=90,hjust=.5,vjust=.5,face="plain"),axis.title.y.right = ggplot2::element_text(colour="grey20",size=30,angle=90,hjust=.5,vjust=.5,face="plain"))


  for(i in 1:dim(rf_data_frame)[1]){

    temprf_data_frame <- data.frame(x=c(0,1),y=c(rf_data_frame$approx_PAF[i],rf_data_frame$approx_PAF[i]))
    p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y), linetype="dashed", size=1.5, col='blue')
  }
  p
}





