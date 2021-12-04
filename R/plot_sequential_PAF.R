#' Produce plots of sequential and average PAF
#'
#' @param SAF_summary An R object produced by running the joint_paf function.
#' @param nrows integer How many rows of plots will be included on the associated figure.
#' @return A plot illustrating average sequential PAF by position and average PAF by risk factor.
plot_sequential <- function(SAF_summary,number_rows){

  riskfactors <- unique(gsub(pattern="(.*)_[0-9]*$",replacement="\\1",x=rownames(SAF_summary)[grep(pattern="^.*_[0-9]*$",x=rownames(SAF_summary),perl=TRUE)]))
  average_PAF <- SAF_summary[grep(pattern=paste0("^Average PAF.*$"),x=rownames(SAF_summary)),3]
  riskfactors <- riskfactors[order(average_PAF,decreasing=TRUE)]
  for(i in 1:length(riskfactors)){

     data_i <- SAF_summary[grep(pattern=paste0(riskfactors[i],"_[0-9]*$"),x=rownames(SAF_summary)),]
     colnames(data_i)[3:5] <- c("value","LB","UB")
     data_i$position <- 1:nrow(data_i)
     data_i$type <- "Sequential"
     data_average <- data_i
     data_average$LB <- SAF_summary[grep(pattern=paste0("Average PAF ", riskfactors[i]),x=rownames(SAF_summary)),4]
     data_average$UB <- SAF_summary[grep(pattern=paste0("Average PAF ", riskfactors[i]),x=rownames(SAF_summary)),5]
     data_average$type <- "Average"
     data_i <- rbind(data_i, data_average)
     p_i <- ggplot2::ggplot(data=data_i, ggplot2::aes(x=position, y=value,  colour="red")) + ggplot2::theme_classic()+ ggplot2::geom_point(,size=4)+ggplot2::geom_ribbon(ggplot2::aes(ymin = LB, ymax = UB, fill=type),alpha=0.2,width= 0.5)+ ggplot2::scale_x_continuous("position",breaks=1:nrow(data_i))+ggplot2::scale_y_continuous("Sequential PAF")+ ggplot2::theme(legend.position = "none")+ ggplot2::annotate(geom="text", x=6, y=0.4, label=riskfactors[i],color="black",size=5)
     eval(parse(text=paste0("p",i,"<- p_i")))
  }
  thetext <- paste0("gridExtra::grid.arrange(p1")
  if(length(riskfactors)>=2) for(i in 2:length(riskfactors)) thetext <- paste0(thetext,",p",i)
  thetext <- paste0(thetext,",nrow=",number_rows,")")
  eval(parse(text=thetext))
}

