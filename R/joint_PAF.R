average_paf_no_CI <- function(data, model_list, parent_list, node_vec,  prev=.09, nperm=NULL, correct_order=3, alpha=0.05,vars=NULL, exact=TRUE){
  #oldw <- getOption("warn")
  #options(warn = -1)

  response_col <- (1:length(colnames(data)))[colnames(data) %in% node_vec[length(node_vec)]]
  if(!c("weights") %in% colnames(data)) data$weights = rep(1, nrow(data))
  if(!is.null(prev)){
    w = prev*as.numeric(data[,response_col]==1) + (1-prev)*as.numeric(data[,response_col]==0)
    data$weights=w
  }
  w <- data$weights
  col_list <- numeric(length(node_vec))
  N <- length(col_list)-1
  sim_disease_current_population <- predict(model_list[[N+1]],type="response")

  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list[indexes]
    N <- length(col_list)-1

  }


   if(exact) correct_order=NULL  # skip if exact calculation
  if(!is.null(correct_order)){


    nperm_new <- factorial(N)/(factorial(N-correct_order))

    repeat_n <- 1

    if(is.null(nperm)){
      nperm <- nperm_new
    }
    if(nperm < nperm_new) nperm <- nperm_new

    if(nperm_new < nperm){

      repeat_n <- floor(nperm/nperm_new)
      nperm <- nperm_new*repeat_n

    }

    perm_mat <- matrix(0,nrow=nperm_new,ncol=N)
    perm_mat[,1:correct_order] <- gtools::permutations(N,correct_order)
    perm_mat_temp <- perm_mat
    if(repeat_n >1){
      for(j in 1:repeat_n){

        perm_mat_temp <- rbind(perm_mat_temp,perm_mat)

      }
    }
    perm_mat <- perm_mat_temp
    rm(perm_mat_temp)
    #print(paste0("doing ", nperm, " permutations"))
  }

  if(exact){

    perm_mat <- matrix(ncol=N)
    for(i in 1:N){
      combos <- gtools::combinations(N,i)
      perm_mat <- rbind(perm_mat,cbind(combos,matrix(0,nrow=nrow(combos),ncol=N-i)))

    }
    perm_mat <- perm_mat[-1,]
      nperm <- nrow(perm_mat)
        theorder <- apply(perm_mat,1,order_fun)
    perm_mat <- perm_mat[order(theorder,decreasing=FALSE),]
  }

  joint_PAF_vec <- numeric(nperm)

  SAF_mat <- matrix(0,nrow=nperm,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nperm,ncol=N)
  order_mat <- matrix(0,nrow=nperm,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nperm,ncol=N)


  for(i in 1:nperm){

if(!exact){
    if(is.null(correct_order)) the_order <- col_list[1:N][sample(1:N,N)]
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      if(correct_order < N) the_order[(correct_order+1):N] <- sample(other_indexes,N-correct_order)
      if(N-correct_order==1) the_order[(correct_order+1):N] <- other_indexes
      the_order <- col_list[1:N][the_order]
    }
    reverse_order <- numeric(N)
    for(j in 1:N) reverse_order[j] <- (1:N)[the_order==col_list[j]]

    current_mat <- data
    current_mat_2 <- data
    SAF <- numeric(N)
    SAF_2 <- numeric(N)
    no_intervention <- sim_disease_current_population


    for(j in 1:length(the_order)){

      current_mat <- sim_outnode(data,the_order[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
      current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
      SAF[j] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
      no_intervention <- current_mat[,col_list[N+1]]

    }
    SAF <- SAF/sum(w*sim_disease_current_population)
    SAF_mat[i,] <- SAF[reverse_order]
    order_mat[i,] <- the_order
    reverse_order_mat[i,] <- reverse_order
    if(i %% 100 == 0){
      flush.console()
      print(i)
    }
}
    if(exact){
      # calculations are for joint PAFs rather than sequential PAFs
      # First check permutation to see if it's the same as previous permutation
      no_intervention <- sim_disease_current_population

      start_again=TRUE
      if(i==1){
        old_perm <- rep(0,N)
        number_rf_new <- sum(perm_mat[i,]!=0)
      }
      if(i > 1){
        old_perm <- perm_mat[i-1,]
        number_rf_new <- sum(perm_mat[i,]!=0)
        number_rf_old <- sum(old_perm!=0)
        if((number_rf_new==number_rf_old+1) && all(old_perm[1:number_rf_old]==perm_mat[i,1:number_rf_old])) start_again=FALSE
      }
      if(start_again==FALSE){

        current_mat <- sim_outnode(data,col_list[1:N][perm_mat[i,number_rf_new]],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
        current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
        joint_PAF_vec[i] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
      }
      if(start_again==TRUE){
        current_mat <- data
        for(j in 1:number_rf_new){

          current_mat <- sim_outnode(data,col_list[1:N][perm_mat[i,j]],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)

        }
        current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
        joint_PAF_vec[i] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
      }

    }

  }

  if(exact){
    joint_PAF_vec <- joint_PAF_vec/sum(w*no_intervention)
    SAF_mat_exact <- matrix(0,nrow=N,ncol=N)
    rownames(SAF_mat_exact) <- paste('riskfactor ',1:N)
    colnames(SAF_mat_exact) <- paste('position ',1:N)
    for(i in 1:N){ # risk factor i
      for(j in 1:N){ # position j

        if(j < N) rows_to_look_at <- (1:nperm)[apply(perm_mat[,1:j,drop=FALSE],1,function(x){any(x==i)}) & perm_mat[,j]>0 & perm_mat[,j+1]==0]
        if(j == N) rows_to_look_at <- (1:nperm)[perm_mat[,N]>0]
        for(k in 1:length(rows_to_look_at)){
          joint_PAF_match_row <- 0
          if(j > 1){
          match_row <- perm_mat[rows_to_look_at[k],]
          match_row <- setdiff(match_row,i)
          match_row <- match_row[1:(j-1)]
          match_row <- (1:nperm)[apply(perm_mat,1,function(x){all(x[1:(j-1)]==match_row)&all(x[j:N]==0)})]
          joint_PAF_match_row <- joint_PAF_vec[match_row]
          }
          SAF_mat_exact[i,j] <- ((k-1)/k)*SAF_mat_exact[i,j]+(joint_PAF_vec[rows_to_look_at[k]]-joint_PAF_match_row)/k
          }
      }
    }


    average_PAF <- apply(SAF_mat_exact,1,function(x){mean(x)})
    SAF_mat_exact <- t(SAF_mat_exact)
    colnames(SAF_mat_exact) <- colnames(data)[col_list][1:N]
    names(average_PAF) <- colnames(data)[col_list][1:N]
    #print(list(SAF_mat=SAF_mat_exact,average_PAF=average_PAF,joint_PAF=joint_PAF_vec[N]))
    #options(warn = oldw)
    thedframe <- data.frame(position=c(rep(paste("elimination position",1:N),N),rep("Average",N),"Joint"),"risk factor"=c(rep(colnames(SAF_mat_exact),times=rep(N,N)),colnames(SAF_mat_exact),""),estimate=c(as.vector(SAF_mat_exact),average_PAF,joint_PAF_vec[N]),check.names=FALSE)
    #print(thedframe)
    avepafs <- structure(thedframe,class="SAF_summary")
      return(avepafs)
  }

  colnames(SAF_mat) <- colnames(data)[col_list][1:N]
  colnames(reverse_order_mat) <- colnames(data)[col_list][1:N]

  average_paf=apply(SAF_mat,2,mean)
  joint_paf <- mean(apply(SAF_mat,1,sum))
  SAF_summary <- matrix(0,nrow=N,ncol=N)

  for(i in 1:N){
    for(j in 1:N){
      SAF_summary[i,j] <- mean(SAF_mat[,j][order_mat[,i]==col_list[j]])
    }
  }
  colnames(SAF_summary) <- names(average_paf)
  rownames(SAF_summary) <- paste("elimination position ", (1:N),sep='')

  ME_SAF_summary <- matrix(0,nrow=N,ncol=N)
  colnames(ME_SAF_summary) <- colnames(SAF_mat)

  for(i in 1:N){
    for(j in 1:N){
      ME_SAF_summary[i,j] <- qt(1-alpha/2, df=sum(order_mat[,i]==col_list[j])-1)*sd(SAF_mat[,j][order_mat[,i]==col_list[j]])/sqrt(sum(order_mat[,i]==col_list[j]))
    }
  }
  temp1 <- reshape2::melt(SAF_summary)
  SAF_summary <- cbind(reshape2::melt(SAF_summary),ME=reshape2::melt(ME_SAF_summary)[,3])

  UB2 <- SAF_summary$value+SAF_summary$ME
  LB2 <- SAF_summary$value-SAF_summary$ME

  SAF_summary$LB <- c(LB2)
  SAF_summary$UB <- c(UB2)
  newdf <- data.frame(Var1=rep("Average",N),Var2=names(average_paf),value=as.numeric(average_paf), ME=qt(1-alpha/2, df=nperm-1)*apply(SAF_mat,2,sd)/sqrt(nperm), LB=as.numeric(average_paf)-qt(1-alpha/2, df=nperm-1)*apply(SAF_mat,2,sd)/sqrt(nperm),UB=as.numeric(average_paf)+qt(1-alpha/2, df=nperm-1)*apply(SAF_mat,2,sd)/sqrt(nperm))
  newdf2 <- data.frame(Var1=c("Joint"),Var2=c(""),value=as.numeric(joint_paf),ME=qt(1-alpha/2, df=nperm-1)*sd(apply(SAF_mat,1,sum))/sqrt(nperm), LB=as.numeric(joint_paf)-qt(1-alpha/2, df=nperm-1)*sd(apply(SAF_mat,1,sum))/sqrt(nperm),UB=as.numeric(joint_paf)+qt(1-alpha/2, df=nperm-1)*sd(apply(SAF_mat,1,sum))/sqrt(nperm))

  SAF_summary <- rbind(SAF_summary, newdf,newdf2)
  rownames(SAF_summary) = NULL
  colnames(SAF_summary) <- c("position", "risk factor", "estimate", "Margin error", "lower bound", "Upper bound")
  #print(SAF_summary)
  SAF_summary <- structure(data.frame(SAF_summary,check.names=FALSE),class="SAF_summary")
  #options(warn = oldw)
  return(SAF_summary)

}

#' Calculation of average and sequential paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of fitted models corresponding for the outcome variables in node_vec, with parents as described in parent_vec.  This list must be in the same order as node_vec and parent_list.  Models can be linear (lm), logistic (glm) or ordinal logistic (polr). Non-linear effects of variables (if necessary) should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec (Note that the variable names should be columns in data)
#' @param node_vec A character vector corresponding to the nodes in the Bayesian network (The variable names should be column names in data).  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param prev numeric.  Prevalence of disease.  Only relevant to set for case control datasets.
#' @param exact logical.  Default TRUE. If TRUE, an efficient calculation is used to calculate average PAF, which enables the average PAF from N! permutations, over all N risk factors to be calculated with only 2^N-1 operations.  If FALSE, permutations are sampled
#' @param nperm  Default NULL Number of random permutations used to calculate average and sequential PAF.  If correct_order is set to an integer value, nperm is reset to an integer multiple of factorial(N)/factorial(N-correct_order) depending on the size of nperm.  If nperm is NULL or less than factorial(N)/factorial(N-correct_order), factorial(N)/factorial(N-correct_order) permutations will be sampled.  If nperm is larger than factorial(N)/factorial(N-correct_order), nperm will be reset to the smallest integer multiple of factorial(N)/factorial(N-correct_order) less than the input value of nperm
#' @param correct_order Default 3.  This enforces stratified sampling of permutations where the first correct_order positions of the sampled permutations are evenly distributed over the integers 1 ... N, N being the number of risk factors of interest, over the sampled permutations.  The other positions are randomly sampled.  This automatically sets the number of simulations when nperm=NULL.  For interest, if N=10 and correct_order=3, nperm is set to factorial(10)/factorial(10-3) = 720.  This special resampling reduces Monte Carlo variation in estimated average and sequential PAFs.
#' @param vars A subset of risk factors for which we want to calculate average, sequential and joint PAF
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendant nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repeated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @param ci_level Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the confidence interval (when ci=TRUE)
#' @param ci_level_ME Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the margin of error for the point estimate (only revelant when ci=FALSE and exact=FALSE)
#' @return A SAF_summary object with average joint and sequential PAF for all risk factors in node_vec (or alternatively a subset of those risk factors if specified in vars).
#' @export
#'
#' @references Ferguson, J., O’Connell, M. and O’Donnell, M., 2020. Revisiting sequential attributable fractions. Archives of Public Health, 78(1), pp.1-9.
#' @references Ferguson, J., Alvarez-Iglesias, A., Newell, J., Hinde, J. and O’Donnell, M., 2018. Estimating average attributable fractions with confidence intervals for cohort and case–control studies. Statistical methods in medical research, 27(4), pp.1141-1152
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' #  Simulated data on occupational and environmental exposure to chronic cough from Eide, 1995
#' # First specify the causal graph, in terms of the parents of each node.  Then put into a list
#' parent_urban.rural <- c()
#' parent_smoking.category <- c("urban.rural")
#' parent_occupational.exposure <- c("urban.rural")
#' parent_y <- c("urban.rural","smoking.category","occupational.exposure")
#' parent_list <- list(parent_urban.rural, parent_smoking.category,
#'  parent_occupational.exposure, parent_y)
#' # also specify nodes of graph, in order from root to leaves
#' node_vec <- c("urban.rural","smoking.category","occupational.exposure", "y")
#' # specify a model list according to parent_list
#' # here we use the auxillary function 'automatic fit'
#' model_list=automatic_fit(data=Hordaland_data, parent_list=parent_list,
#'  node_vec=node_vec, prev=.09)
#' # By default the function works by stratified simulation of permutations and
#' # subsequent simulation of the incremental interventions on the distribution of risk
#' # factors.  The permutations are stratified so each factor appears equally often in
#' # the first correct_order positions.  correct_order has a default of 2.
#' out <- average_paf(data=Hordaland_data, model_list=model_list, parent_list=parent_list,
#'  node_vec=node_vec, prev=.09, nperm=10,vars = c("urban.rural",
#'  "occupational.exposure"),ci=FALSE)
#'  print(out)
#'
#' \donttest{
#' # More complicated example (slower to run)
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol",
#' "stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet",
#' "smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol","stress",
#' "lipids","waist_hip_ratio","high_blood_pressure",
#' "early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,
#' parent_alcohol,parent_stress,parent_high_blood_pressure,
#' parent_lipids,parent_waist_hip_ratio,parent_early_stage_heart_disease,
#' parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress",
#' "high_blood_pressure","lipids","waist_hip_ratio","early_stage_heart_disease",
#' "diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list,
#'  node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)",
#'   spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' out <- average_paf(data=stroke_reduced, model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.0035,
#' vars = c("high_blood_pressure","smoking","stress","exercise","alcohol",
#' "diabetes","early_stage_heart_disease"),ci=TRUE,boot_rep=10)
#' print(out)
#' plot(out,max_PAF=0.5,min_PAF=-0.1,number_rows=3)
#' # plot sequential and average PAFs by risk factor
#' # similar calculation, but now sampling permutations (stratified, so
#' # that each risk factor will appear equally often in the first correct_order positions)
#' out <- average_paf(data=stroke_reduced, model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.0035, exact=FALSE,
#'  correct_order=2, vars = c("high_blood_pressure","smoking","stress",
#'  "exercise","alcohol","diabetes","early_stage_heart_disease"),ci=TRUE,
#'  boot_rep=10)
#'  print(out)
#'  plot(out,max_PAF=0.5,min_PAF=-0.1,number_rows=3)
#' }
average_paf <- function(data, model_list, parent_list, node_vec, prev=.09, exact=TRUE, nperm=NULL, correct_order=2, vars=NULL,ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95){

  if(!node_order(parent_list=parent_list,node_vec=node_vec)){
    stop("ancestors must be specified before descendants in node_vec")
  }
  if(!is.null(vars) & !all(vars %in% node_vec)){
    stop("Not all requested variables are in node_vec.  Check spelling")
  }
  if(!is.null(correct_order) && is.null(vars)) correct_order <- min(correct_order,length(node_vec))
  if(!is.null(correct_order) && !is.null(vars)) correct_order <- min(correct_order,length(vars))
  if(is.null(correct_order)&&is.null(nperm)){

    stop("please specify either correct_order and nperm")

  }
  ## how many risk factors are under scrutiny?
  col_list <- numeric(length(node_vec))
  N <- length(col_list)-1
  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list[indexes]
    N <- length(col_list)-1

  }

  if(!ci){
    res <- average_paf_no_CI(data=data, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nperm=nperm, correct_order=correct_order, alpha=1-ci_level_ME,vars=vars,exact=exact)
    return(res)
  }
  #oldw <- getOption("warn")
  #options(warn = -1)
  nc <- options()$boot.ncpus
  cl <- parallel::makeCluster(nc)
  parallel::clusterExport(cl, c("ns"))
  res <- boot::boot(data=data,statistic=average_paf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nperm=nperm, correct_order=correct_order, vars=vars, exact=exact,cl=cl)
  parallel::stopCluster(cl)
  if(is.null(vars)) vars <- node_vec[1:(length(node_vec)-1)]
  res <- extract_ci(res=res,model_type='glm',t_vector=c(paste0(rep(node_vec[node_vec %in% vars],times=rep(length(vars),length(vars))),'_',rep(1:length(vars),length(vars))),paste0("Average PAF ", node_vec[node_vec %in% vars]),'JointPAF'),ci_level=ci_level,ci_type=ci_type,continuous=TRUE)
  res <- cbind(position=c(rep(paste("elimination position",1:N),N),rep("Average",N),"Joint"),'risk factor'=rownames(res),res)
  rownames(res) <- NULL
  res$`risk factor` <- gsub(pattern="(.*)_[0-9]",replacement="\\1",x=res$`risk factor`,perl=TRUE)
  res$`risk factor` <- gsub(pattern="Average PAF (.*)",replacement="\\1",x=res$`risk factor`,perl=TRUE)
  res$`risk factor` <- gsub(pattern="Joint",replacement="",x=res$`risk factor`,perl=TRUE)
  #print(res)
  res <- structure(res,class="SAF_summary")
  #options(warn = oldw)
  return(res)

}

#' Print out SAF_summary object
#'
#' @param x A SAF_summary object.  This is a special dataframe that is created by running the function average_PAF.
#' @param ... Other arguments to be passed to print
#'
#' @return No return value.  Prints the SAF_summary object to the console.
#' @export
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' #  Simulated data on occupational and environmental exposure to chronic cough from Eide, 1995
#' # First specify the causal graph, in terms of the parents of each node.  Then put into a list
#' parent_urban.rural <- c()
#' parent_smoking.category <- c("urban.rural")
#' parent_occupational.exposure <- c("urban.rural")
#' parent_y <- c("urban.rural","smoking.category","occupational.exposure")
#'parent_list <- list(parent_urban.rural, parent_smoking.category,
#'parent_occupational.exposure, parent_y)
# # also specify nodes of graph, in order from root to leaves
#' node_vec <- c("urban.rural","smoking.category","occupational.exposure", "y")
# # specify a model list according to parent_list
# # here we use the auxillary function 'automatic fit'
#' model_list=automatic_fit(data=Hordaland_data, parent_list=parent_list,
#'  node_vec=node_vec, prev=.09)
#' out <- average_paf(data=Hordaland_data, model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.09, nperm=10,
#' vars = c("urban.rural","occupational.exposure"),ci=FALSE)
print.SAF_summary <- function(x,...){

  data_frame <- structure(as.list(x),class="data.frame", row.names=attr(x,"row.names"))
  print(data_frame)

}

average_paf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09, nperm=100, correct_order=3, vars=NULL, exact=TRUE){

  ################################

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
      model_text <- paste0("survival::clogit(",model_text,",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    if(model_type=="coxph"){

      model_text <- as.character(model$call)
      model_text <- paste0("survival::coxph(",model_text[2],",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type== "glm"){
      #browser()
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
      if(with_weights==TRUE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=weights)")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type == "lm"){
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(with_weights==TRUE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data,weights=weights)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type == "polr"){
      model_text <- as.character(model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    model
  }



  sim_outnode <- function(data,col_num, current_mat, parent_list, col_list,model_list){

    if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(data[,col_num])[1]
    if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

    colname <- colnames(current_mat)[col_num]

    for(i in 1:(length(parent_list)-1)){
      if(colname %in% parent_list[[i]]){
        if(length(table(current_mat[,col_list[[i]]] ))==1) next

        if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
        if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim(col_list[i],current_mat,model_list[[i]],SN=TRUE)
      }
    }
    current_mat
  }



  do_sim <- function(colnum,current_mat, model,SN=FALSE){
    ## polr
    if(names(model)[2]=='zeta'){

      probs <- predict(model,newdata=current_mat,type="probs")
      mynames <- colnames(probs)
      return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
    }
    # glm
    if(length(grep("glm",model$call))>0){

      probs <- predict(model,newdata=current_mat,type="response")
      if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
      return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
    }
    # regression
    if(length(grep("lm",model$call))>0){

      pred <- predict(model,newdata=current_mat,type="response")
      resids <- model$residuals
      if(SN){

        return(pred+resids)

      }

      #browser()
      #return(pred + sample(summary(model)$residuals,length(resids),replace=TRUE))
      #return(pred + rnorm(length(resids),mean=0,sd=.1*sd(resids)))
      return(pred + sample(resids,length(resids),replace=TRUE, prob=model$weights/sum(model$weights)))
    }
  }
  ##################################


  data <- data[ind,]
  n_data <- nrow(data)
  response_col <- (1:length(colnames(data)))[colnames(data) %in% node_vec[length(node_vec)]]
  if(!c("weights") %in% colnames(data)) data$weights = rep(1, nrow(data))
  if(!is.null(prev)){
    w = prev*as.numeric(data[,response_col]==1) + (1-prev)*as.numeric(data[,response_col]==0)
    data$weights=w
  }
  w <- data$weights
  #  if(!all(ind==1:n_data)) browser()
  if(!all(ind==1:n_data)) for(i in 1:length(model_list)) model_list[[i]] <- refit(model=model_list[[i]],data=data)


   col_list <- numeric(length(node_vec))
  N <- length(col_list)-1
  sim_disease_current_population <- predict(model_list[[N+1]],type="response")

  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list[indexes]
    N <- length(col_list)-1

  }



    if(exact) correct_order=NULL  # skip if exact calculation
  if(!is.null(correct_order)){

    nperm_new <- factorial(N)/(factorial(N-correct_order))

    repeat_n <- 1

    if(is.null(nperm)){
      nperm <- nperm_new
    }
    if(nperm < nperm_new) nperm <- nperm_new

    if(nperm_new < nperm){

      repeat_n <- floor(nperm/nperm_new)
      nperm <- nperm_new*repeat_n

    }

    perm_mat <- matrix(0,nrow=nperm_new,ncol=N)
    perm_mat[,1:correct_order] <- gtools::permutations(N,correct_order)
    perm_mat_temp <- perm_mat
    if(repeat_n >1){
    for(j in 1:repeat_n){

      perm_mat_temp <- rbind(perm_mat_temp,perm_mat)

    }
    }
    perm_mat <- perm_mat_temp
    rm(perm_mat_temp)

  }

  order_fun <- function(x){

    N <- length(x)
    sum <- 0
    for(i in 1:N){
      sum <- sum + x[i]*(N+1)^(N-i)
    }
    return(sum)
  }
  if(exact){

    perm_mat <- matrix(ncol=N)
    for(i in 1:N){
      combos <- gtools::combinations(N,i)
      perm_mat <- rbind(perm_mat,cbind(combos,matrix(0,nrow=nrow(combos),ncol=N-i)))

    }
    perm_mat <- perm_mat[-1,]
     nperm <- nrow(perm_mat)
       theorder <- apply(perm_mat,1,order_fun)
    perm_mat <- perm_mat[order(theorder,decreasing=FALSE),]
    }


  SAF_mat <- matrix(0,nrow=nperm,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nperm,ncol=N)
  order_mat <- matrix(0,nrow=nperm,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nperm,ncol=N)
  joint_PAF_vec <- numeric(nperm) # only used when exact
   for(i in 1:nperm){

  if(!exact){
    if(is.null(correct_order)) the_order <- col_list[1:N][sample(1:N,N)]
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      if(correct_order < N) the_order[(correct_order+1):N] <- sample(other_indexes,N-correct_order)
      if(N-correct_order==1) the_order[(correct_order+1):N] <- other_indexes
      the_order <- col_list[1:N][the_order]
    }

    reverse_order <- numeric(N)
    for(j in 1:N) reverse_order[j] <- (1:N)[the_order==col_list[j]]

    current_mat <- data
    current_mat_2 <- data
    SAF <- numeric(N)
    SAF_2 <- numeric(N)
    no_intervention <- sim_disease_current_population


    for(j in 1:length(the_order)){

      current_mat <- sim_outnode(data,the_order[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
      current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
      SAF[j] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
      no_intervention <- current_mat[,col_list[N+1]]

    }
    SAF <- SAF/sum(w*sim_disease_current_population)
    SAF_mat[i,] <- SAF[reverse_order]
    order_mat[i,] <- the_order
    reverse_order_mat[i,] <- reverse_order
  }
   if(exact){
     # calculations are for joint PAFs rather than sequential PAFs
     # First check permutation to see if it's the same as previous permutation
     no_intervention <- sim_disease_current_population

     start_again=TRUE
     if(i==1){
       old_perm <- rep(0,N)
       number_rf_new <- sum(perm_mat[i,]!=0)
     }
     if(i > 1){
       old_perm <- perm_mat[i-1,]
        number_rf_new <- sum(perm_mat[i,]!=0)
        number_rf_old <- sum(old_perm!=0)
        if((number_rf_new==number_rf_old+1) && all(old_perm[1:number_rf_old]==perm_mat[i,1:number_rf_old])) start_again=FALSE
     }
     if(start_again==FALSE){

       current_mat <- sim_outnode(data,col_list[1:N][perm_mat[i,number_rf_new]],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
       current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
       joint_PAF_vec[i] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
     }
     if(start_again==TRUE){
       current_mat <- data
       for(j in 1:number_rf_new){

         current_mat <- sim_outnode(data,col_list[1:N][perm_mat[i,j]],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)

       }
       current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
       joint_PAF_vec[i] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
     }

}


  }


  if(exact){
    joint_PAF_vec <- joint_PAF_vec/sum(w*no_intervention)
    SAF_mat_exact <- matrix(0,nrow=N,ncol=N)
    rownames(SAF_mat_exact) <- paste('riskfactor ',1:N)
    colnames(SAF_mat_exact) <- paste('position ',1:N)
    for(i in 1:N){ # risk factor i
      for(j in 1:N){ # position j

        if(j < N) rows_to_look_at <- (1:nperm)[apply(perm_mat[,1:j,drop=FALSE],1,function(x){any(x==i)}) & perm_mat[,j]>0 & perm_mat[,j+1]==0]
        if(j == N) rows_to_look_at <- (1:nperm)[perm_mat[,N]>0]
        for(k in 1:length(rows_to_look_at)){
          joint_PAF_match_row <- 0
          if(j > 1){
            match_row <- perm_mat[rows_to_look_at[k],]
            match_row <- setdiff(match_row,i)
            match_row <- match_row[1:(j-1)]
            match_row <- (1:nperm)[apply(perm_mat,1,function(x){all(x[1:(j-1)]==match_row)&all(x[j:N]==0)})]
            joint_PAF_match_row <- joint_PAF_vec[match_row]
          }
          SAF_mat_exact[i,j] <- ((k-1)/k)*SAF_mat_exact[i,j]+(joint_PAF_vec[rows_to_look_at[k]]-joint_PAF_match_row)/k
        }
      }
    }

    average_PAF <- apply(SAF_mat_exact,1,function(x){mean(x)})
    SAF_mat_exact <- t(SAF_mat_exact)
    colnames(SAF_mat_exact) <- colnames(data)[col_list][1:N]
    names(average_PAF) <- colnames(data)[col_list][1:N]
    return(c(SAF_mat=as.numeric(SAF_mat_exact),average_PAF=average_PAF,joint_PAF=joint_PAF_vec[N]))
  }

  colnames(SAF_mat) <- colnames(data)[col_list][1:N]
  colnames(reverse_order_mat) <- colnames(data)[col_list][1:N]

  average_paf=apply(SAF_mat,2,mean)
  joint_paf <- mean(apply(SAF_mat,1,sum))
  SAF_summary <- matrix(0,nrow=N,ncol=N)

  for(i in 1:N){
    for(j in 1:N){
      SAF_summary[i,j] <- mean(SAF_mat[,j][order_mat[,i]==col_list[j]])
    }
  }
  return_vec <- c(as.numeric(SAF_summary),average_paf,joint_paf)
  return(return_vec)

}



sim_outnode <- function(data,col_num, current_mat, parent_list, col_list,model_list){

  if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(data[,col_num])[1]
  if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

  colname <- colnames(current_mat)[col_num]

  for(i in 1:(length(parent_list)-1)){
    if(colname %in% parent_list[[i]]){
      if(length(table(current_mat[,col_list[[i]]] ))==1) next

      if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
      if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim(col_list[i],current_mat,model_list[[i]],SN=TRUE)
    }
  }
  current_mat
}



do_sim <- function(colnum,current_mat, model,SN=FALSE){
  ## polr
  if(names(model)[2]=='zeta'){

    probs <- predict(model,newdata=current_mat,type="probs")
    mynames <- colnames(probs)
    return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
  }
  # glm
  if(length(grep("glm",model$call))>0){

    probs <- predict(model,newdata=current_mat,type="response")
    if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
    return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
  }
  # regression
  if(length(grep("lm",model$call))>0){

    pred <- predict(model,newdata=current_mat,type="response")
    resids <- model$residuals
    if(SN){

      return(pred+resids)

    }

    #browser()
    #return(pred + sample(summary(model)$residuals,length(resids),replace=TRUE))
    #return(pred + rnorm(length(resids),mean=0,sd=.1*sd(resids)))
    return(pred + sample(resids,length(resids),replace=TRUE, prob=model$weights/sum(model$weights)))
  }
}

make_formula <- function(parents,outcome_node,common='',spline_nodes=c(),df_spline_nodes=3){
   if(length(parents)==0) return(paste(outcome_node,"~ 1"))
  if(common!="") result <- paste(outcome_node,"~",common,"+ ",parents[1])
  if(common=="") result <- paste(outcome_node,"~",parents[1])
  if(length(parents)>=2){

    for(i in 2:length(parents)){

      if(parents[i] %in% spline_nodes) result <- paste(result,"+ns(",parents[i],",df=",df_spline_nodes,')',sep='')
      if(!parents[i] %in% spline_nodes) result <- paste(result,"+ ",parents[i],sep='')

    }
  }
  result
}
#' Automatic fitting models for Bayesian network.
#'
#' Main effects models are fit by default.  For continuous variables, lm is used, for binary (numeric 0/1 variables), glm is used and for factor valued variables polr is used.  For factors, ensure that the factor levels are ordered by increasing levels of risk.  If interactions are required for certain models, it is advisable to populate the elements of model_list separately.
#'
#' @param data Data frame. A data frame containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param prev  Prevalence of disease.  Set to NULL for cohort or cross sectional studies
#' @param common character text for part of the model formula that doesn't involve any variable in node_vec.  Useful for specifying confounders involved in all models automatically
#' @param spline_nodes  Vector of continuous variable names that are fit as splines (when involved as parents).  Natural splines are used.
#' @param df_spline_nodes How many degrees of freedom for each spline (Default 3).  At the moment, this can not be specified separately for differing variables.
#' @return A list of fitted models corresponding to node_vec and parent_vec.
#' @export
#'
#' @examples
#' # More complicated example (slower to run)
#' library(splines)
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet",
#' "smoking","alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet",
#' "smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol","stress",
#' "lipids","waist_hip_ratio","high_blood_pressure","early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,
#' parent_alcohol,parent_stress,parent_high_blood_pressure,
#' parent_lipids,parent_waist_hip_ratio,parent_early_stage_heart_disease,
#' parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure",
#' "lipids","waist_hip_ratio","early_stage_heart_disease",
#' "diabetes","case")
#' \donttest{
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list,
#' node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+
#' sex*ns(age,df=5)", spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' }
automatic_fit <- function(data, parent_list, node_vec, prev=.09,common='',spline_nodes=c(),df_spline_nodes=3){
  #oldw <- getOption("warn")
  #options(warn = -1)

model_list=list()
outcome_name <- node_vec[length(node_vec)]
outcome_bin <- data[,colnames(data) %in% outcome_name]
if(!is.null(prev)){
data$weights=1
data$weights = prev*as.numeric(outcome_bin==1) + (1-prev)*as.numeric(outcome_bin==0)
}
if(!c("weights") %in% colnames(data)) data$weights <- rep(1,nrow(data))
for(i in 1:length(node_vec)){
  column <- (1:length(colnames(data)))[colnames(data) %in% node_vec[i]]
  formula_text <- make_formula(parent_list[[i]],node_vec[i],common=common,spline_nodes=spline_nodes,df_spline_nodes=df_spline_nodes)
  thesplit <- ""
  while(length(grep(pattern='^.*ns\\(.*$',x=formula_text))>0){
    formula_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=formula_text)
    stuff <- strsplit(formula_text,split="splines::ns(",fixed=TRUE)
    formula_text <- stuff[[1]][1]
    thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
  }
  formula_text <- paste0(formula_text,thesplit)

  y <- data[,column]

  if(i < length(node_vec)){
        if(length(table(y))==2){
      theform <- paste("glm(",formula_text,",data=data,family='binomial',weights=weights)",sep='')
    }
    if(length(table(y))>2 & is.factor(y)){
      theform <- paste("MASS::polr(",formula_text,",data=data,weights=weights)",sep='')
    }
    if(length(table(y))>2 & is.numeric(y)){
      theform <- paste("lm(",formula_text,",data=data,weights=weights)",sep='')
    }
  }
  if(i==length(node_vec)) theform <- paste("glm(",formula_text,",data=data,family='binomial',weights=weights)",sep='')

  to_execute <- paste("model_list[[i]] <-", theform,sep='')
  eval(parse(text=to_execute))
}
#options(warn = oldw)

model_list
}

node_order <- function(parent_list, node_vec){

  L <- length(node_vec)
  for(i in 1:(L-1)){

    putative_ancestors <- unique(unlist(parent_list[1:i]))
    if(any(node_vec[(i+1):L] %in% putative_ancestors)) return(FALSE)

  }
  return(TRUE)
}

order_fun <- function(x){

  N <- length(x)
  sum <- 0
  for(i in 1:N){
    sum <- sum + x[i]*(N+1)^(N-i)
  }
  return(sum)
}


##################################  the same functions as above are replicated here - but only return joint_PAF

#' Calculation of joint paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of fitted models corresponding for the outcome variables in node_vec, with parents as described in parent_vec.  This list must be in the same order as node_vec and parent_list. Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.  Linear (lm), logistic (glm) and ordinal logistic (polr) models are permitted
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param prev prevalence of the disease (default is NULL)
#' @param vars A subset of risk factors for which we want to calculate joint PAF
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendant nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repeated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available method
#' @param ci_level Numeric.  Confidence level.  Default 0.95
#' @param nsim Numeric.  Number of independent simulations of the dataset.  Default of 1.
#' @return A numeric estimate of the joint PAF for all risk factors (if ci=FALSE), or a data frame giving joint PAF and confidence intervals (if ci=TRUE)
#' @export
#'
#' @references Ferguson, J., O’Connell, M. and O’Donnell, M., 2020. Revisiting sequential attributable fractions. Archives of Public Health, 78(1), pp.1-9.
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#' #  Simulated data on occupational and environmental exposure to
#' # chronic cough from Eide, 1995
#' # First specify the causal graph, in terms of the parents of each node.
#' # Then put into a list.
#' parent_urban.rural <- c()
#' parent_smoking.category <- c("urban.rural")
#' parent_occupational.exposure <- c("urban.rural")
#' parent_y <- c("urban.rural","smoking.category","occupational.exposure")
#' parent_list <- list(parent_urban.rural, parent_smoking.category,
#'  parent_occupational.exposure, parent_y)
#' # also specify nodes of graph, in order from root to leaves
#' node_vec <- c("urban.rural","smoking.category","occupational.exposure", "y")
#' # specify a model list according to parent_list
#' # here we use the auxillary function 'automatic fit'
#' model_list=automatic_fit(data=Hordaland_data, parent_list=parent_list,
#' node_vec=node_vec, prev=.09)
#' joint_paf(data=Hordaland_data, model_list=model_list, parent_list=parent_list,
#'  node_vec=node_vec, prev=.09, vars = c("urban.rural",
#'  "occupational.exposure"),ci=FALSE)
#' \donttest{
#' # More complicated example (slower to run)
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet",
#' "smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure","early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,parent_alcohol,
#' parent_stress,parent_high_blood_pressure,parent_lipids,parent_waist_hip_ratio,
#' parent_early_stage_heart_disease,parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure",
#' "lipids","waist_hip_ratio","early_stage_heart_disease","diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list,
#' node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)",
#'  spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' jointpaf <- joint_paf(data=stroke_reduced, model_list=model_list,
#' parent_list=parent_list, node_vec=node_vec, prev=.0035,
#' vars = c("high_blood_pressure","smoking","stress","exercise","alcohol",
#' "diabetes","early_stage_heart_disease"),ci=TRUE,boot_rep=10)
#' }
joint_paf <- function(data, model_list, parent_list, node_vec, prev=NULL, vars=NULL,ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95,nsim=1){
  if(!node_order(parent_list=parent_list,node_vec=node_vec)){
    stop("ancestors must be specified before descendants in node_vec")
  }
  if(!is.null(vars) & !all(vars %in% node_vec)){
    stop("Not all requested variables are in node_vec.  Check spelling")
  }
if(!ci) return(joint_paf_inner(data=data,ind=1:nrow(data), model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev,vars=vars,nsim=nsim))
  nc <- options()$boot.ncpus
  cl <- parallel::makeCluster(nc)
  parallel::clusterExport(cl, c("ns"))

  res <- boot::boot(data=data,statistic=joint_paf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, vars=vars,nsim=nsim,cl=cl)
  parallel::stopCluster(cl)
  stuff <- extract_ci(res=res,model_type='glm',ci_level=ci_level,ci_type=ci_type,continuous=TRUE,t_vector=c("joint PAF"))
  return(stuff)

}


joint_paf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09,vars=NULL,nsim=1){

   ################################

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
      model_text <- paste0("survival::clogit(",model_text,",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    if(model_type=="coxph"){

      model_text <- as.character(model$call)
      model_text <- paste0("survival::coxph(",model_text[2],",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type== "glm"){
      #browser()
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
      if(with_weights==TRUE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=weights)")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)

    }


      model <- eval(parse(text=model_text))
    }

    if(model_type == "lm"){
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(with_weights==TRUE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data,weights=weights)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type == "polr"){
      model_text <- as.character(model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    model
  }


  sim_outnode <- function(data,col_num, current_mat, parent_list, col_list,model_list){

    if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(data[,col_num])[1]
    if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

    colname <- colnames(current_mat)[col_num]

    for(i in 1:(length(parent_list)-1)){
      if(colname %in% parent_list[[i]]){
        if(length(table(current_mat[,col_list[[i]]] ))==1) next

        if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
        if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim(col_list[i],current_mat,model_list[[i]],SN=TRUE)
      }
    }
    current_mat
  }



  do_sim <- function(colnum,current_mat, model,SN=FALSE){
    ## polr
    if(names(model)[2]=='zeta'){

      probs <- predict(model,newdata=current_mat,type="probs")
      mynames <- colnames(probs)
      return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
    }
    # glm
    if(length(grep("glm",model$call))>0){

      probs <- predict(model,newdata=current_mat,type="response")
      if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
      return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
    }
    # regression
    if(length(grep("lm",model$call))>0){

      pred <- predict(model,newdata=current_mat,type="response")
      resids <- model$residuals
      if(SN){

        return(pred+resids)

      }

      #browser()
      #return(pred + sample(summary(model)$residuals,length(resids),replace=TRUE))
      #return(pred + rnorm(length(resids),mean=0,sd=.1*sd(resids)))
      return(pred + sample(resids,length(resids),replace=TRUE, prob=model$weights/sum(model$weights)))
    }
  }
  ##################################


  data <- data[ind,]
  n_data <- nrow(data)
  response_col <- (1:length(colnames(data)))[colnames(data) %in% node_vec[length(node_vec)]]
  if(!c("weights") %in% colnames(data)) data$weights = rep(1, nrow(data))
  if(!is.null(prev)){
    w = prev*as.numeric(data[,response_col]==1) + (1-prev)*as.numeric(data[,response_col]==0)
    data$weights=w
  }
  w <- data$weights
  #  if(!all(ind==1:n_data)) browser()
  if(!all(ind==1:n_data)) for(i in 1:length(model_list)) model_list[[i]] <- refit(model=model_list[[i]],data=data)


   sim_disease_current_population <- predict(model_list[[length(node_vec)]],type="response")

  out_vector <- numeric(nsim)
  for(k in 1:nsim){
    col_list <- numeric(length(node_vec))
    N <- length(col_list)-1
  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list_orig[indexes]
    N <- length(col_list)-1

  }
current_mat <- data
      for(j in 1:N){

        current_mat <- sim_outnode(data,col_list[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
        current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")

      }
  out_vector[k] <- as.numeric((sum(w*sim_disease_current_population)-sum(w*current_mat[,col_list[N+1]]))/sum(w*sim_disease_current_population))
}
      return(jointPAF=mean(out_vector))

}

#' Calculation of Sequential paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of fitted model objects corresponding for the outcome variables in node_vec, with parents as described in parent_vec. Linear (lm), logistic (glm) and ordinal (polr) objects are allowed. This list must be in the same order as node_vec and parent_list.  Non-linear effects should be specified via ns(x, df=y), where ns is the natural spline function from the splines library.
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param prev prevalence of the disease (default is NULL)
#' @param vars A character vector of riskfactors.  Sequential PAF is calculated for the risk factor specified in the last position of the vector, conditional on the other risk factors
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendant nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repeated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE).  Note that at least 50 replicates are recommended to achieve stable estimates of standard error.  In the examples below, values of boot_rep less than 50 are sometimes used to limit run time.
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @param ci_level Numeric.  Confidence level.  Default 0.95
#' @param nsim Numeric.  Number of independent simulations of the dataset.  Default of 1.
#' @return A numeric estimate of sequential PAF (if ci=FALSE), or else a data frame giving estimates and confidence limits of sequential PAF (if ci=TRUE)
#' @export
#'
#' @references Ferguson, J., O’Connell, M. and O’Donnell, M., 2020. Revisiting sequential attributable fractions. Archives of Public Health, 78(1), pp.1-9.
#'
#' @examples
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=2)
#' # The above could be set to the number of available cores on the machine
#'
#' # Simulated data on occupational and environmental exposure to
#' # chronic cough from Eide, 1995
#' # First specify the causal graph, in terms of the parents of each node.
#' # Then put into a list.
#' parent_urban.rural <- c()
#' parent_smoking.category <- c("urban.rural")
#' parent_occupational.exposure <- c("urban.rural")
#' parent_y <- c("urban.rural","smoking.category","occupational.exposure")
#' parent_list <- list(parent_urban.rural, parent_smoking.category,
#'  parent_occupational.exposure, parent_y)
#' # also specify nodes of graph, in order from root to leaves
#' node_vec <- c("urban.rural","smoking.category","occupational.exposure", "y")
#' # specify a model list according to parent_list
#' # here we use the auxillary function 'automatic fit'
#' model_list=automatic_fit(data=Hordaland_data, parent_list=parent_list,
#' node_vec=node_vec, prev=.09)
#' # sequential paf for occupational exposure conditional on elimination of urban.rural
#' seq_paf(data=Hordaland_data, model_list=model_list, parent_list=parent_list,
#'  node_vec=node_vec, prev=.09, vars = c("urban.rural",
#'  "occupational.exposure"),ci=FALSE)
#' \donttest{
#' # More complicated example (slower to run)
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking","alcohol",
#' "stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking",
#' "alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet",
#' "smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol",
#' "stress","lipids","waist_hip_ratio","high_blood_pressure",
#' "early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,parent_alcohol,
#' parent_stress,parent_high_blood_pressure,parent_lipids,parent_waist_hip_ratio,
#' parent_early_stage_heart_disease,parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure",
#' "lipids","waist_hip_ratio","early_stage_heart_disease","diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list,
#' node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)",
#'  spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' # calculate sequential PAF for stress, conditional on smoking
#' # and blood pressure being eliminated from the population
#' seqpaf <- seq_paf(data=stroke_reduced, model_list=model_list, parent_list=
#' parent_list, node_vec=node_vec, prev=.0035, vars = c("high_blood_pressure",
#' "smoking","stress"),ci=TRUE,boot_rep=10)
#' }
seq_paf <- function(data, model_list, parent_list, node_vec, prev=NULL, vars=NULL,ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95,nsim=1){
  if(!node_order(parent_list=parent_list,node_vec=node_vec)){
    stop("ancestors must be specified before descendants in node_vec")
  }
  if(!is.null(vars) & !all(vars %in% node_vec)){
    stop("Not all requested variables are in node_vec.  Check spelling")
  }
  if(!is.null(vars) & length(vars)<2){
    stop("Enter at least 2 risk factors.  SAF is calculated for the last risk factor conditional on the others in list")
  }
  if(!ci) return(seq_paf_inner(data=data,ind=1:nrow(data), model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev,vars=vars,nsim=nsim))
  nc <- options()$boot.ncpus
  cl <- parallel::makeCluster(nc)
  parallel::clusterExport(cl, c("ns"))
  res <- boot::boot(data=data,statistic=seq_paf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, vars=vars,nsim=nsim,cl=cl)
  parallel::stopCluster(cl)
  stuff <- extract_ci(res=res,model_type='glm',ci_level=ci_level,ci_type=ci_type,continuous=TRUE,t_vector=c("sequential PAF"))
  return(stuff)

}


seq_paf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09,vars=NULL,nsim=1){



  ################################


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
      model_text <- paste0("survival::clogit(",model_text,",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    if(model_type=="coxph"){

      model_text <- as.character(model$call)
      model_text <- paste0("survival::coxph(",model_text[2],",data=data)")
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type== "glm"){
      #browser()
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
      if(with_weights==TRUE && length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=weights)")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type == "lm"){
      model_text <- as.character(model$call)
      if(with_weights==FALSE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(with_weights==TRUE && length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data,weights=weights)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }

    if(model_type == "polr"){
      model_text <- as.character(model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      model_text <- model_text_u
      if(length(grep(pattern='^(.*)splines::ns\\((.*)$',x=model_text))==0){
        thesplit <- ""
        while(length(grep(pattern='^.*ns\\(.*$',x=model_text))>0){
          model_text <- gsub(pattern='^(.*)ns\\((.*)$',replacement='\\1splines::ns\\(\\2',x=model_text)
          stuff <- strsplit(model_text,split="splines::ns(",fixed=TRUE)
          model_text <- stuff[[1]][1]
          thesplit <- paste0("splines::ns(",stuff[[1]][2],thesplit)
        }
        model_text <- paste0(model_text,thesplit)
      }
      model <- eval(parse(text=model_text))
    }
    model
  }


  sim_outnode <- function(data,col_num, current_mat, parent_list, col_list,model_list){

    if(is.factor(current_mat[,col_num])) current_mat[,col_num] <- levels(data[,col_num])[1]
    if(is.numeric(current_mat[,col_num])) current_mat[,col_num] <- 0

    colname <- colnames(current_mat)[col_num]

    for(i in 1:(length(parent_list)-1)){
      if(colname %in% parent_list[[i]]){
        if(length(table(current_mat[,col_list[[i]]] ))==1) next

        if(is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- factor(do_sim(col_list[i],current_mat,model_list[[i]]),levels=levels(current_mat[,col_list[i]]))
        if(!is.factor(current_mat[,col_list[i]])) current_mat[,col_list[i]] <- do_sim(col_list[i],current_mat,model_list[[i]],SN=TRUE)
      }
    }
    current_mat
  }



  do_sim <- function(colnum,current_mat, model,SN=FALSE){
    ## polr
    if(names(model)[2]=='zeta'){

      probs <- predict(model,newdata=current_mat,type="probs")
      mynames <- colnames(probs)
      return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
    }
    # glm
    if(length(grep("glm",model$call))>0){

      probs <- predict(model,newdata=current_mat,type="response")
      if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
      return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
    }
    # regression
    if(length(grep("lm",model$call))>0){

      pred <- predict(model,newdata=current_mat,type="response")
      resids <- model$residuals
      if(SN){

        return(pred+resids)

      }

      #browser()
      #return(pred + sample(summary(model)$residuals,length(resids),replace=TRUE))
      #return(pred + rnorm(length(resids),mean=0,sd=.1*sd(resids)))
      return(pred + sample(resids,length(resids),replace=TRUE, prob=model$weights/sum(model$weights)))
    }
  }
  ##################################


  data <- data[ind,]
  n_data <- nrow(data)
  response_col <- (1:length(colnames(data)))[colnames(data) %in% node_vec[length(node_vec)]]
  if(!c("weights") %in% colnames(data)) data$weights = rep(1, nrow(data))
  if(!is.null(prev)){
    w = prev*as.numeric(data[,response_col]==1) + (1-prev)*as.numeric(data[,response_col]==0)
    data$weights=w
  }
  w <- data$weights
  #  if(!all(ind==1:n_data)) browser()
  if(!all(ind==1:n_data)) for(i in 1:length(model_list)) model_list[[i]] <- refit(model=model_list[[i]],data=data)


  sim_disease_current_population <- predict(model_list[[length(node_vec)]],type="response")

  out_vector <- numeric(nsim)
   for(k in 1:nsim){
    col_list <- numeric(length(node_vec))
    N <- length(col_list)-1
    for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
    col_list_orig <- col_list
    if(!is.null(vars)){
      #browser()
      indexes <- numeric(length(vars))
      for(i in 1:length(vars))  indexes[i] <- (1:(N+1))[node_vec %in% vars[i]]
      indexes <- c(indexes,N+1)
      col_list <- col_list_orig[indexes]
      N <- length(col_list)-1

    }
    current_mat <- data
    for(j in 1:(N-1)){

      current_mat <- sim_outnode(data,col_list[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
      current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")

    }
    PAF_S <- as.numeric((sum(w*sim_disease_current_population)-sum(w*current_mat[,col_list[N+1]]))/sum(w*sim_disease_current_population))
      current_mat <- sim_outnode(data,col_list[N],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
    current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
    PAF_S_J <- as.numeric((sum(w*sim_disease_current_population)-sum(w*current_mat[,col_list[N+1]]))/sum(w*sim_disease_current_population))
    out_vector[k] <- PAF_S_J-PAF_S
  }
  return(jointPAF=mean(out_vector))

}


