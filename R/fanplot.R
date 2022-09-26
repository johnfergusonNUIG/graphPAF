#' Create a rf.data.frame object
#'
#' Create a rf.data.frame object for risk factors, prevalence and risk ratios.  This will be used in fan plots and nomograms (by simply sending the rf.dat.frame object to plot)
#'
#' @param rf_names A character vector of risk factor names
#' @param rf_prev A numeric vector specifying prevalence of risk factor in disease controls (estimates of population prevalence can also be used if the disease is rare)
#' @param risk A numeric vector of relative risks or Odds ratios for disease corresponding to each risk factor (if log=FALSE).  Log-relative risks or log-odds ratios can  be alternatively specified (if log=TRUE)
#' @param log default TRUE. Set to TRUE if relative risks/odds ratios are specified on log-scale
#' @return A rf.data.frame object
#' @export
#'
#' @references Ferguson, J., O’Leary, N., Maturo, F., Yusuf, S. and O’Donnell, M., 2019. Graphical comparisons of relative disease burden across multiple risk factors. BMC medical research methodology, 19(1), pp.1-9.
#'
#' @examples
#' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA','Diet',
#' 'WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),
#' rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),
#' risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)
#' # fanplot
#' plot(rfs,type="f")
#' # nomogram
#' plot(rfs,type="n")
#' # reverse nomogram
#' # plot(rfs,type="rn")
rf_summary <- function(rf_names, rf_prev, risk, log=FALSE){
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
#' @param x A rf.data.frame object
#' @param type A character representing the type of plot.  "f" for a fan_plot, "n" for a PAF nomogram and "rn" for a reverse PAF nomogram
#' @param rf_prevmarks Axis marks for risk factor prevalence (only used for type="n" and type = "rn") Default c(0.02, 0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
#' @param ormarks Axis marks for odds ratios (only used for type="n" and type = "rn") Default c(1.05,1.1,1.4,1.7,2.0,3.0)
#' @param type A character representing the type of plot.  "f" for a fan_plot, "n" for a PAF nomogram and "rn" for a reverse PAF nomogram.  See Ferguson et al.. "Graphical comparisons of relative disease burden across multiple risk factors." BMC medical research methodology 19, no. 1 (2019): 1-9 for more details
#' @param ...  Other arguments that can be passed to the plotting routine
#' @return fanplot or PAF nomogram (each is a ggplotting object)
#' @export
#'
#' @references Ferguson, J., O’Leary, N., Maturo, F., Yusuf, S. and O’Donnell, M., 2019. Graphical comparisons of relative disease burden across multiple risk factors. BMC medical research methodology, 19(1), pp.1-9.
#'
#' @examples
#' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA',
#' 'Diet','WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),
#' rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),
#' risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)
#' # fanplot
#' plot(rfs,type="f")
#' # nomogram
#' plot(rfs,type="n")
#' # reverse nomogram
#' plot(rfs,type="rn")
plot.rf.data.frame <- function(x,type="f", rf_prevmarks= c(0.02, 0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9),ormarks = c(1.05,1.1,1.4,1.7,2.0,3.0), ...)
{
    rf_data_frame <- x
    if(!inherits(rfs,"rf.data.frame")){

    stop("Create a valid rf.data.frame object before running function")

  }
  rf_data_frame <- structure(as.list(rf_data_frame),class="data.frame", row.names=attr(rf_data_frame,"row.names"))
  if(type=="f"){

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

  }
 if(type=="n"){
   ggplot2::theme_set(ggplot2::theme_classic())
   a <- max(max(-1*log(rf_data_frame$rf_prev)+0.1*abs(log(rf_data_frame$rf_prev))),max(log(rf_data_frame$approx_PAF)))
   b <-  min(min(-1*log(rf_data_frame$rf_prev)-0.1*abs(log(rf_data_frame$rf_prev))),min(log(rf_data_frame$approx_PAF)))

   s <- (b*max(-log(rf_prevmarks))-a*min(-log(rf_prevmarks)))/(b-a)
   c <- a/(max(-log(rf_prevmarks))-s)

   newylimits <- 1.08*c(min(c(c*(-log(rf_prevmarks)-s),c*(log(rf_prevmarks)+s))),max(c(c*(-log(rf_prevmarks)-s),c*(log(rf_prevmarks)+s))))



   # Plot
   p <- ggplot2::ggplot(rf_data_frame) + ggplot2::geom_segment(ggplot2::aes(x=0.5, xend=2.5, y=c*(-log(rf_prev)-s), yend=c*(log(approx_PAF)+s), col=rf_data_frame$rf_names), size=.75) +
     # color of lines
     ggplot2::xlim(0, 3) + ggplot2::ylim(newylimits)  # X and Y axis limits

   p <- p + ggplot2::labs(col = "")


   for(i in 1:length(rf_prevmarks)){

     temprf_data_frame <- data.frame(x=c(0.4,0.5),y=rep(c*(-log(rf_prevmarks)[i]-s),2))
     p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y), linetype="dashed")
     temprf_data_frame <- data.frame(x=c(2.5,2.6),y=rep(c*(log(rf_prevmarks)[i]+s),2))
     p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y), linetype="dashed")
     temprf_data_frame <- data.frame(x=c(1.5,1.55),y=rep(.5*c*log(log(ormarks))[i],2))
     p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y), linetype="dashed")
   }

   ## add axes

   temprf_data_frame <- data.frame(x=c(1.5,1.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y))

   temprf_data_frame <- data.frame(x=c(0.5,0.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y))

   temprf_data_frame <- data.frame(x=c(2.5,2.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=temprf_data_frame,ggplot2::aes(x=x,y=y))
   ## add tickmarks


   # Add texts

   p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(0.35, 9),y=c*(-log(rf_prevmarks)-s), size=6, fontface=2)
   p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(2.65, 9),y=c*(log(c(0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))+s), size=6,fontface=2)
   p <- p + ggplot2::annotate(geom="text",label=c("1.05","1.1","1.4","1.7","2.0","3.0"),x=rep(1.6, 6),y=.5*c*log(log(c(1.05,1.1,1.4,1.7,2.0,3))), size=6, fontface=2)

   p <- p + ggrepel::geom_label_repel(ggplot2::aes(label=rf_data_frame$rf_names,x=rep(0.5, nrow(rf_data_frame)),y=c*(-log(rf_prev)-s)), size=6, data=rf_data_frame)
   #p <- p + geom_text(aes(label=1:nrow(rf_data_frame),x=rep(2.5, NROW(rf_data_frame)),y=log(approx_PAF)), size=4, data=rf_data_frame)
   p <- p + ggrepel::geom_label_repel(ggplot2::aes(label=rf_data_frame$rf_names,x=rep(2.5, nrow(rf_data_frame)),y=c*(log(approx_PAF)+s)), size=6, data=rf_data_frame)

   p <- p + ggplot2::geom_text(label="approximate PAF", x=2.5, y=min(min(log(rf_data_frame$rf_prev)-0.1*abs(log(rf_data_frame$rf_prev))),min(log(rf_data_frame$approx_PAF)-0.1*abs(log(rf_data_frame$approx_PAF)))), hjust=0.5, size=5)  # title

   p <- p + ggplot2::geom_text(label="Odds Ratio", x=1.5, y=min(min(log(rf_data_frame$rf_prev)-0.1*abs(log(rf_data_frame$rf_prev))),min(log(rf_data_frame$approx_PAF)-0.1*abs(log(rf_data_frame$approx_PAF)))), hjust=0.5, size=5)  # title

   p <- p + ggplot2::geom_text(label="Prevalence in Controls", x=0.5, y=min(min(log(rf_data_frame$rf_prev)-0.1*abs(log(rf_data_frame$rf_prev))),min(log(rf_data_frame$approx_PAF)-0.1*abs(log(rf_data_frame$approx_PAF)))), hjust=0.5, size=5)  # title

   # Minimal theme
   p <- p + ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),axis.line = ggplot2::element_blank(),axis.ticks = ggplot2::element_blank(),axis.text.x = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank(),axis.title.x = ggplot2::element_blank(),axis.title.y= ggplot2::element_blank(),panel.border = ggplot2::element_blank(), ggplot2::element_blank(),plot.margin = ggplot2::unit(c(1,2,1,2), "cm"),legend.text=ggplot2::element_text(size=16))


 }
 if(type=="rn"){

   a <- max(max(-1*log(rf_data_frame$log_riskratio)+0.02*abs(log(rf_data_frame$log_riskratio))),max(log(rf_data_frame$approx_PAF)+0.02*abs(log(rf_data_frame$approx_PAF))))
   b <-  min(min(-1*log(rf_data_frame$log_riskratio)-0.02*abs(log(rf_data_frame$log_riskratio))),min(log(rf_data_frame$approx_PAF)-0.02*abs(log(rf_data_frame$approx_PAF))))

   a1 <- max(max(-1*log(rf_data_frame$log_riskratio)+0.1*abs(log(rf_data_frame$log_riskratio))),max(log(rf_data_frame$approx_PAF)+0.1*abs(log(rf_data_frame$approx_PAF))))
   b1 <-  min(min(-1*log(rf_data_frame$log_riskratio)-0.1*abs(log(rf_data_frame$log_riskratio))),min(log(rf_data_frame$approx_PAF)-0.1*abs(log(rf_data_frame$approx_PAF))))

   newylimits=c(b1,a1)



   s1 <- (a-b)/(max(log(rf_prevmarks))-min(log(rf_prevmarks)))
   s2 <- (a-b)/(max(-log(log(ormarks)))-min(-log(log(ormarks))))

   # scale
   s <- min(s1,s2)

   prev_add <- (a+b)/2  - s*((max(log(rf_prevmarks))+min(log(rf_prevmarks))))/2
   or_add <- (a+b)/2 - s*(max(-log(log(ormarks)))+min(-log(log(ormarks))))/2

   p <- ggplot2::ggplot(rf_data_frame) + ggplot2::geom_segment(ggplot2::aes(x=0.5, xend=2.5, y=s*(-log(rf_data_frame$log_riskratio))+or_add, yend=s*(log(rf_data_frame$approx_PAF))+prev_add, col=rf_data_frame$rf_names), size=.75) +
     # color of lines
     ggplot2::xlim(0, 3)+ggplot2::ylim(newylimits)  # X and Y axis limits

   p <- p + ggplot2::labs(col = "")


   for(i in 1:length(rf_prevmarks)){

     tempdf <- data.frame(x=c(0.4,0.5),y=rep(s*(-log(log(ormarks)))[i]+or_add),2)
     p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y), linetype="dashed")
     tempdf <- data.frame(x=c(2.5,2.6),y=rep(s*(log(rf_prevmarks)[i])+prev_add,2))
     p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y), linetype="dashed")
     tempdf <- data.frame(x=c(1.5,1.55),y=rep(.5*s*log(rf_prevmarks)[i]+or_add*.5+prev_add*.5,2))
     p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y), linetype="dashed")
   }

   ## add axes

   tempdf <- data.frame(x=c(1.5,1.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y))

   tempdf <- data.frame(x=c(0.5,0.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y))

   tempdf <- data.frame(x=c(2.5,2.5),y=newylimits)
   p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y))

   # Add texts ....

   p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(1.6, length(rf_prevmarks)),y=.5*(s*(log(rf_prevmarks)))+.5*(prev_add)+.5*or_add, size=6, fontface=2)
   p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(2.65, length(rf_prevmarks)),y=s*(log(rf_prevmarks))+prev_add, size=6,fontface=2)
   p <- p + ggplot2::annotate(geom="text",label=ormarks,x=rep(.35, length(ormarks)),y=s*(-log(log(ormarks)))+or_add, size=6, fontface=2)

   p <- p + ggrepel::geom_label_repel(ggplot2::aes(label=rf_data_frame$rf_names,x=rep(0.5, nrow(rf_data_frame)),y=s*(-log(rf_data_frame$log_riskratio))+or_add), size=6, data=rf_data_frame)
   #p <- p + geom_text(aes(label=order,x=rep(2.5, NROW(df)),y=log(rf_data_frame$approx_paf)), size=6, data=df)
   p <- p + ggrepel::geom_label_repel(ggplot2::aes(label=rf_data_frame$rf_names,x=rep(2.5, nrow(rf_data_frame)),y=s*(log(rf_data_frame$approx_PAF))+prev_add), size=6, data=rf_data_frame)

   p <- p + ggplot2::geom_text(label="approximate PAF", x=2.5, y=b-.08*(a-b), hjust=0.5, size=5)  # title

   p <- p + ggplot2::geom_text(label="Prevalence in Controls", x=1.5, y=b-.08*(a-b), hjust=0.5, size=5)  # title

   p <- p + ggplot2::geom_text(label="Odds Ratio", x=0.5, y=b-.08*(a-b), hjust=0.5, size=5)  # title

   # Minify theme
   p <- p + ggplot2::theme(panel.background = ggplot2::element_blank(),
                           panel.grid = ggplot2::element_blank(),axis.line = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y= ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           plot.margin = ggplot2::unit(c(1,2,1,2), "cm"),legend.text=ggplot2::element_text(size=16))

 }
  p
}





