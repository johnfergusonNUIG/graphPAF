#' Produce plots of sequential and average PAF
#'
#' @param x An SAF_summary R object produced by running the average_paf function.
#' @param number_rows integer How many rows of plots will be included on the associated figure.
#' @param max_PAF upper limit of y axis on PAF plots (default = 0.4)
#' @param min_PAF lower limit of y axis on PAF plots (default = 0)
#' @param ... Other global arguments inherited by that might be passed to the ggplot routine
#' @return A ggplot2 plotting object illustrating average sequential PAF by position and average PAF by risk factor.
#'
#' @references Ferguson, J., O’Connell, M. and O’Donnell, M., 2020. Revisiting sequential attributable fractions. Archives of Public Health, 78(1), pp.1-9.
#' Ferguson, J., Alvarez-Iglesias, A., Newell, J., Hinde, J. and O’Donnell, M., 2018. Estimating average attributable fractions with confidence intervals for cohort and case–control studies. Statistical methods in medical research, 27(4), pp.1141-1152
#'
#' @export
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' #  Simulated data on occupational and environmental exposure to
#' # chronic cough from Eide, 1995
#' # First specify the causal graph, in terms of the parents of each node.  Then put into a list
#' parent_urban.rural <- c()
#' parent_smoking.category <- c("urban.rural")
#' parent_occupational.exposure <- c("urban.rural")
#' parent_y <- c("urban.rural","smoking.category","occupational.exposure")
#' parent_list <- list(parent_urban.rural, parent_smoking.category,
#' parent_occupational.exposure, parent_y)
#' # also specify nodes of graph, in order from root to leaves
#' node_vec <- c("urban.rural","smoking.category","occupational.exposure", "y")
#' model_list=automatic_fit(Hordaland_data,
#' parent_list=parent_list, node_vec=node_vec, prev=.09)
#' out <- average_paf(data=model_list[[length(model_list)]]$data,
#'  model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.09, nperm=10,
#' vars = c("urban.rural","occupational.exposure"),ci=FALSE)
#' plot(out)
#' \donttest{
#' # plot with confidence intervals for average and sequential PAF
#' # (This is probably more useful for more than 2 risk factors).
#' # Separate axes for each risk factor so confidence intervals can be clearly displayed
#' out <- average_paf(data=model_list[[length(model_list)]]$data,
#'  model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.09, nperm=10,
#' vars = c("urban.rural","occupational.exposure"),ci=TRUE,boot_rep=8)
#' plot(out)
#' # Here we plot, with margin of error of point estimate when 50 permutations are used
#' out <- average_paf(data=model_list[[length(model_list)]]$data,
#'  model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.09, nperm=50,
#' vars = c("urban.rural","occupational.exposure"),ci=FALSE,exact=FALSE)
#' plot(out)
#' }
plot.SAF_summary <- function(x,number_rows=3, max_PAF=0.4,min_PAF=0,...){
  #defaultW <- getOption("warn")
  #options(warn = -1)
  SAF_summary <- x
  SAF_summary <- structure(as.list(SAF_summary),class="data.frame", row.names=attr(SAF_summary,"row.names"))

  riskfactors <- unique(SAF_summary$`risk factor`)
  average_PAF <- SAF_summary[grep(pattern=paste0("^Average.*$"),SAF_summary$position),3]
  riskfactors <- riskfactors[order(average_PAF,decreasing=TRUE)]
  N <- ncol(SAF_summary)
  if("Margin error"%in%colnames(SAF_summary)){
    SAF_summary <- SAF_summary[,-grep("Margin error",colnames(SAF_summary))]
    N <- N-1
  }
  if(N>3){
  for(i in 1:length(riskfactors)){

     data_i <- SAF_summary[intersect(grep(pattern=paste("elimination position"),x=SAF_summary$position),grep(pattern=riskfactors[i],x=SAF_summary$`risk factor`)),]
       colnames(data_i)[(N-2):N] <- c("value","LB","UB")
     data_i$position <- 1:nrow(data_i)
     data_i$type <- "Sequential"
     data_average <- data_i
     data_average$value <- SAF_summary[intersect(grep(pattern=paste("Average"),x=SAF_summary$position),grep(pattern=riskfactors[i],x=SAF_summary$`risk factor`)),N-2]
     data_average$LB <- SAF_summary[intersect(grep(pattern=paste("Average"),x=SAF_summary$position),grep(pattern=riskfactors[i],x=SAF_summary$`risk factor`)),N-1]
     data_average$UB <- SAF_summary[intersect(grep(pattern=paste("Average"),x=SAF_summary$position),grep(pattern=riskfactors[i],x=SAF_summary$`risk factor`)),N]
     data_average$type <- "Average"
     data_i <- rbind(data_i, data_average)
     p_i <- ggplot2::ggplot(data=data_i, ggplot2::aes(x=position, y=value,  colour=type)) + ggplot2::theme_classic()+ ggplot2::geom_point(,size=4)+ggplot2::geom_ribbon(ggplot2::aes(ymin = LB, ymax = UB, fill=type),alpha=0.2,width= 0.5)+ ggplot2::scale_x_continuous("position",breaks=1:nrow(data_i))+ggplot2::scale_y_continuous("Sequential PAF",limits=c(min_PAF,max_PAF))+ ggplot2::theme(legend.position = "none")+ ggplot2::annotate(geom="text", x = quantile(data_i$position,.5), y = max_PAF, label=riskfactors[i],color="black",size=5)
     eval(parse(text=paste0("p",i,"<- p_i")))
  }
  thetext <- paste0("gridExtra::grid.arrange(p1")
  if(length(riskfactors)>=2) for(i in 2:length(riskfactors)) thetext <- paste0(thetext,",p",i)
  thetext <- paste0(thetext,",nrow=",number_rows,")")
  eval(parse(text=thetext))
  }
  if(N==3){

          data_elim <- SAF_summary[grep(pattern=paste("elimination position"),x=SAF_summary$position),]
      colnames(data_elim)[N] <- c("value")
      data_elim$position <- as.numeric(gsub(pattern="elimination position (.*)",replacement="\\1",data_elim$position))
      p_i <- ggplot2::ggplot(data=data_elim, ggplot2::aes(x=position,y=value,color=`risk factor`))+ ggplot2::theme_classic()+ggplot2::geom_line(ggplot2::aes(color=`risk factor`))+ ggplot2::scale_x_continuous("position",breaks=1:nrow(data_elim))+ggplot2::scale_y_continuous("Sequential PAF",limits=c(min_PAF,max_PAF))+ggplot2::labs(title="Sequential PAF")+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      p_i

  }
  #options(warn = defaultW)
}

