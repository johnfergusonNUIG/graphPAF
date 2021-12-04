

#' Create a PAF nomogram of a rf.data.frame object
#'
#' @param rf_data_frame
#'
#' @return
#' @export
#'
#' @examples
#'
#' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA','Diet','WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)
#' nomogram(rfs)

nomogram <- function(rf_data_frame){

   if(!class(rf_data_frame)=='rf.data.frame'){

    stop("Create a valid rf.data.frame object before running function")

  }
  rf_data_frame <- structure(as.list(rf_data_frame),class="data.frame", row.names=attr(rf_data_frame,"row.names"))


ggplot2::theme_set(ggplot2::theme_classic())
a <- max(max(-1*log(rf_data_frame$rf_prev)+0.1*abs(log(rf_data_frame$rf_prev))),max(log(rf_data_frame$approx_PAF)))
b <-  min(min(-1*log(rf_data_frame$rf_prev)-0.1*abs(log(rf_data_frame$rf_prev))),min(log(rf_data_frame$approx_PAF)))

rf_prevmarks <- c(0.02, 0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
ormarks <- c(1.05,1.1,1.4,1.7,2.0,3.0)

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
p
}


#' Nomogram with Odds ratio/Risk Ratio on Left hand aix
#'
#' @param rf_data_frame
#' @param prevmarks
#' @param ormarks
#'
#' @return
#' @export
#'
#' @examples
#'
#' #' rfs <- rf_summary(rf_names=c('Hypertension','Inactivity','ApoB/ApoA','Diet','WHR','Smoking','Cardiac causes','Alcohol','Global Stress','Diabetes'),rf_prev=c(.474,.837,.669,.67,.67,.224,.049,.277,.144,.129),risk=c(1.093,0.501,0.428,0.378,0.294,0.513,1.156,0.186,0.301,0.148),log=TRUE)
#' reverse_nomogram(rfs)


reverse_nomogram <- function(rf_data_frame,prevmarks = c(0.02, 0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9),ormarks=c(1.05,1.1,1.2,1.3,1.5,1.7,2.0,2.5,3.0)){


  if(!class(rf_data_frame)=='rf.data.frame'){

    stop("Create a valid rf.data.frame object before running function")

  }

  rf_data_frame <- structure(as.list(rf_data_frame),class="data.frame", row.names=attr(rf_data_frame,"row.names"))

  rf_prevmarks <- c(0.02, 0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
  ormarks <- c(1.05,1.1,1.4,1.7,2.0,3.0)

    a <- max(max(-1*log(rf_data_frame$log_riskratio)+0.02*abs(log(rf_data_frame$log_riskratio))),max(log(rf_data_frame$approx_PAF)+0.02*abs(log(rf_data_frame$approx_PAF))))
    b <-  min(min(-1*log(rf_data_frame$log_riskratio)-0.02*abs(log(rf_data_frame$log_riskratio))),min(log(rf_data_frame$approx_PAF)-0.02*abs(log(rf_data_frame$approx_PAF))))

    a1 <- max(max(-1*log(rf_data_frame$log_riskratio)+0.1*abs(log(rf_data_frame$log_riskratio))),max(log(rf_data_frame$approx_PAF)+0.1*abs(log(rf_data_frame$approx_PAF))))
    b1 <-  min(min(-1*log(rf_data_frame$log_riskratio)-0.1*abs(log(rf_data_frame$log_riskratio))),min(log(rf_data_frame$approx_PAF)-0.1*abs(log(rf_data_frame$approx_PAF))))

    newylimits=c(b1,a1)



    s1 <- (a-b)/(max(log(prevmarks))-min(log(prevmarks)))
    s2 <- (a-b)/(max(-log(log(ormarks)))-min(-log(log(ormarks))))

    # scale
    s <- min(s1,s2)

    prev_add <- (a+b)/2  - s*((max(log(prevmarks))+min(log(prevmarks))))/2
    or_add <- (a+b)/2 - s*(max(-log(log(ormarks)))+min(-log(log(ormarks))))/2

    p <- ggplot2::ggplot(rf_data_frame) + ggplot2::geom_segment(ggplot2::aes(x=0.5, xend=2.5, y=s*(-log(rf_data_frame$log_riskratio))+or_add, yend=s*(log(rf_data_frame$approx_PAF))+prev_add, col=rf_data_frame$rf_names), size=.75) +
      # color of lines
      ggplot2::xlim(0, 3)+ggplot2::ylim(newylimits)  # X and Y axis limits

    p <- p + ggplot2::labs(col = "")


    for(i in 1:length(rf_prevmarks)){

      tempdf <- data.frame(x=c(0.4,0.5),y=rep(s*(-log(log(ormarks)))[i]+or_add),2)
      p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y), linetype="dashed")
      tempdf <- data.frame(x=c(2.5,2.6),y=rep(s*(log(prevmarks)[i])+prev_add,2))
      p <- p + ggplot2::geom_line(data=tempdf,ggplot2::aes(x=x,y=y), linetype="dashed")
      tempdf <- data.frame(x=c(1.5,1.55),y=rep(.5*s*log(prevmarks)[i]+or_add*.5+prev_add*.5,2))
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

    p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(1.6, length(prevmarks)),y=.5*(s*(log(prevmarks)))+.5*(prev_add)+.5*or_add, size=6, fontface=2)
    p <- p + ggplot2::annotate(geom="text",label=c("2%","5%","10%","20%","30%","40%","50%","70%","90%"),x=rep(2.65, length(prevmarks)),y=s*(log(prevmarks))+prev_add, size=6,fontface=2)
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
    p

}
