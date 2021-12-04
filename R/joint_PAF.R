joint_paf_no_CI <- function(data, model_list, parent_list, node_vec,  prev=.09, nsim=NULL, correct_order=3, alpha=0.05,vars=NULL){
  response_col <- (1:length(colnames(data)))[colnames(data) %in% node_vec[length(node_vec)]]
  if(!c("weights") %in% colnames(data)) data$weights = rep(1, nrow(data))
  if(!is.null(prev)){
    w = prev*as.numeric(data[,response_col]==1) + (1-prev)*as.numeric(data[,response_col]==0)
    data$weights=w
  }
  w <- data$weights
  col_list <- numeric(length(node_vec))
  N <- length(col_list)-1

  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list[indexes]

  }

  sim_disease_current_population <- predict(model_list[[N+1]],type="response")

  N <- length(col_list)-1

  if(!is.null(correct_order)){


    nsim_new <- factorial(N)/(factorial(N-correct_order))

    repeat_n <- 1

    if(is.null(nsim)){
      nsim <- nsim_new
    }
    if(nsim < nsim_new) nsim <- nsim_new

    if(nsim_new < nsim){

      repeat_n <- floor(nsim/nsim_new)
      nsim <- nsim_new*repeat_n

    }

    perm_mat <- matrix(0,nrow=nsim_new,ncol=N)
    perm_mat[,1:correct_order] <- gtools::permutations(N,correct_order)
    perm_mat_temp <- perm_mat
    if(repeat_n >1){
      for(j in 1:repeat_n){

        perm_mat_temp <- rbind(perm_mat_temp,perm_mat)

      }
    }
    perm_mat <- perm_mat_temp
    rm(perm_mat_temp)
    print(paste0("doing ", nsim, " permutations"))
  }

  SAF_mat <- matrix(0,nrow=nsim,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nsim,ncol=N)
  order_mat <- matrix(0,nrow=nsim,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nsim,ncol=N)


  for(i in 1:nsim){

    if(is.null(correct_order)) the_order <- col_list[1:N][sample(1:N,N)]
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      if(correct_order < N) the_order[(correct_order+1):N] <- sample(other_indexes,N-correct_order)
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
  newdf <- data.frame(Var1=rep("Average",N),Var2=names(average_paf),value=as.numeric(average_paf), ME=qt(1-alpha/2, df=nsim-1)*apply(SAF_mat,2,sd)/sqrt(nsim), LB=as.numeric(average_paf)-qt(1-alpha/2, df=nsim-1)*apply(SAF_mat,2,sd)/sqrt(nsim),UB=as.numeric(average_paf)+qt(1-alpha/2, df=nsim-1)*apply(SAF_mat,2,sd)/sqrt(nsim))
  newdf2 <- data.frame(Var1=c("Joint"),Var2=c(""),value=as.numeric(joint_paf),ME=qt(1-alpha/2, df=nsim-1)*sd(apply(SAF_mat,1,sum))/sqrt(nsim), LB=as.numeric(joint_paf)-qt(1-alpha/2, df=nsim-1)*sd(apply(SAF_mat,1,sum))/sqrt(nsim),UB=as.numeric(joint_paf)+qt(1-alpha/2, df=nsim-1)*sd(apply(SAF_mat,1,sum))/sqrt(nsim))

  SAF_summary <- rbind(SAF_summary, newdf,newdf2)
  rownames(SAF_summary) = NULL
  colnames(SAF_summary) <- c("position", "risk factor", "estimate", "Margin error", "lower bound", "Upper bound")
  return(SAF_summary)

}

#' Calculation of joint, average and sequential paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of models corresponding for the outcome variables in node_vec, with parents as described in parent_vec.  This list must be in the same order as node_vec and parent_list
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param nsim  Default NULL Number of random permutations used to calculate average and sequential PAF.  If correct_order is set to an integer value, nsim is reset to the largest integer multiple of correct_order that is less than the number of permutations implied by correct_order.
#' @param correct_order Default 3.  This enforces stratified sampling of permutations where the first correct_order positions of the sampled permutations are evenly distributed over the integers 1 ... n, n being the number of risk factors of interest, over the sampled permutations.  The other positions are randomly sampled.  This automatically sets the number of simulations.  For interest, if n=10 and correct_order=3, nsim is set to factorial(n)/factorial(n-correct_order).  This special resampling reduces Monte Carlo variation in estimated average and sequential PAFs.
#' @params vars A subset of risk factors for which we want to calculate average, sequential and joint PAF
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendent nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE)
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @param ci_level Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the confidence interval (when ci=TRUE)
#' @param ci_level_ME Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the margin of error for the point estimate
#' @return A dataframe with average joint and sequential PAF for all risk factors in node_vec (or alternatively a subset of those risk factors if specified in vars).
#' @export

joint_paf <- function(data, model_list, parent_list, node_vec, prev=.09, nsim=NULL, correct_order=2, vars=NULL,ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95){
  if(!node_order(parent_list=parent_list,node_vec=node_vec)){
    stop("ancestors must be specified before descendants in node_vec")
  }
  if(!is.null(vars) & !all(vars %in% node_vec)){
    stop("Not all requested variables are in node_vec.  Check spelling")
  }
  if(!is.null(correct_order) && is.null(vars)) correct_order <- min(correct_order,length(node_vec))
  if(!is.null(correct_order) && !is.null(vars)) correct_order <- min(correct_order,length(vars))
  if(is.null(correct_order)&&is.null(nsim)){

    stop("please specify either correct_order and nsim")

  }
  if(!ci) return(joint_paf_no_CI(data=data, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nsim=nsim, correct_order=correct_order, alpha=1-ci_level_ME,vars=vars))
  res <- boot::boot(data=data,statistic=joint_paf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nsim=nsim, correct_order=correct_order, vars=vars)
  if(is.null(vars)) vars <- node_vec[1:(length(node_vec)-1)]

      return(extract_ci(res=res,model_type='glm',t_vector=c(paste0(rep(node_vec[node_vec %in% vars],times=rep(length(vars),length(vars))),'_',rep(1:length(vars),length(vars))),paste0("Average PAF ", node_vec[node_vec %in% vars]),'JointPAF'),ci_level=ci_level,ci_type=ci_type,continuous=TRUE))

}

joint_paf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09, nsim=100, correct_order=3, vars=NULL){

  library(splines)
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
  for(i in 1:(N+1)) col_list[i] <- (1:ncol(data))[colnames(data)==node_vec[i]]
  col_list_orig <- col_list
  if(!is.null(vars)){
    #browser()
    indexes <- c((1:(N+1))[node_vec %in% vars],N+1)
    col_list <- col_list[indexes]

  }

  sim_disease_current_population <- predict(model_list[[length(node_vec)]],type="response")

  N <- length(col_list)-1
  if(!is.null(correct_order)){

    nsim_new <- factorial(N)/(factorial(N-correct_order))

    repeat_n <- 1

    if(is.null(nsim)){
      nsim <- nsim_new
    }
    if(nsim < nsim_new) nsim <- nsim_new

    if(nsim_new < nsim){

      repeat_n <- floor(nsim/nsim_new)
      nsim <- nsim_new*repeat_n

    }

    perm_mat <- matrix(0,nrow=nsim_new,ncol=N)
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
  SAF_mat <- matrix(0,nrow=nsim,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nsim,ncol=N)
  order_mat <- matrix(0,nrow=nsim,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nsim,ncol=N)


  for(i in 1:nsim){

    if(is.null(correct_order)) the_order <- col_list[1:N][sample(1:N,N)]
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      if(correct_order < N) the_order[(correct_order+1):N] <- sample(other_indexes,N-correct_order)
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

      if(parents[i] %in% spline_nodes) result <- paste(result,"+ ns(",parents[i],",df=",df_spline_nodes,')',sep='')
      if(!parents[i] %in% spline_nodes) result <- paste(result,"+ ",parents[i],sep='')

    }
  }
  result
}
automatic_fit <- function(data, parent_list, node_vec, prev=.09,common='',spline_nodes=c(),df_spline_nodes=3){

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





