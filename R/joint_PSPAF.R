average_pspaf_no_CI <- function(data, model_list, parent_list, node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = NULL, alpha=0.05, vars = NULL, exact = TRUE, response_model, mediator_models, riskfactor, refval, calculation_method = "D"){

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

  if(exact){

    perm_mat <- matrix(ncol=N)
    for(i in 1:N){
      combos <- gtools::combinations(N,i)
      perm_mat <- rbind(perm_mat,cbind(combos,matrix(0,nrow=nrow(combos),ncol=N-i)))

    }
    perm_mat <- perm_mat[-1,]
      nsim <- nrow(perm_mat)
        theorder <- apply(perm_mat,1,order_fun)
    perm_mat <- perm_mat[order(theorder,decreasing=FALSE),]
  }

  joint_PAF_vec <- numeric(nsim)

  SAF_mat <- matrix(0,nrow=nsim,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nsim,ncol=N)
  order_mat <- matrix(0,nrow=nsim,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nsim,ncol=N)

######################################################################################
######
## Begininng of material to MOVE above for loop
######
######################################################################################
  Num_forboot <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- NULL
  if(grepl("^glm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "glm"
  if(grepl("^lm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "lm"
  if(grepl("^.*polr$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "polr"

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    if(grepl("^glm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "glm"
    if(grepl("^lm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "lm"
    if(grepl("^.*polr$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "polr"
  }

###################################################################################################
mediatorNames <- as.data.frame( lapply(mediator_models,function(x){as.character(formula(x))[2]}))

pathOrder <- cbind(riskfactor,mediatorNames, as.character( formula(response_model) )[2])

pathspecific_col_list <- numeric(length(pathOrder))

  N <- length(pathspecific_col_list)-1

  for(i in 1:(N+1)) pathspecific_col_list[i] <- (1:ncol(data))[colnames(data)==pathOrder[1,i]]

  pathspecific_col_list_orig <- pathspecific_col_list

  riskfactor_col <- (1:length(colnames(data)))[colnames(data) %in% riskfactor]
######################################################################################
  ######
  ## End of material to MOVE above for loop
  ######
######################################################################################

  for(i in 1:nsim){

if(!exact){

    if(is.null(correct_order)){
      the_order <- pathspecific_col_list[1:N][sample(1:N,N)]
      the_order_colNums1toN <- numeric(N)
      for(j in 1:N) the_order_colNums1toN[j] <- (1:N)[pathspecific_col_list_orig==the_order[j] ]
    }
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      if(correct_order < N) the_order[(correct_order+1):N] <- as.numeric( sample(as.character(other_indexes),N-correct_order) )
      the_order_colNums1toN <- the_order
      the_order <- col_list[1:N][the_order]
    }
    reverse_order <- numeric(N)
    for(j in 1:N) reverse_order[j] <- (1:N)[the_order==pathspecific_col_list[j]]

    current_mat <- data
    current_mat_2 <- data
    SAF <- numeric(N)
    SAF_2 <- numeric(N)
    no_intervention <- sim_disease_current_population

    #   ###########################################################
    #   ###########################################################
    #   ## Previous Sequential PAF code
    # for(j in 1:length(the_order)){
    #   # current_mat <- sim_outnode(data,the_order[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
    #   # current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
    #   # SAF[j] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
    #   # no_intervention <- current_mat[,col_list[N+1]]
    # }
    #   ###########################################################
    #   ###########################################################
###################################################################################################
###################################################################################################
### CODE ADD IN FOR Permutations  SEQUENTIAL AND JOINT PSPAF
###################################################################################################
###################################################################################################
###############################
###############################

                  # initialise previous_SAF at zero and same size as output from impact_fraction()
                  previous_SAF <- impact_fraction(model=response_model, data=data, new_data=data,calculation_method=calculation_method, prev=prev,ci=FALSE)
                  for(j in 1:length(the_order)){

                    if( the_order[j]  == riskfactor_col ){

                          current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                    }else{
                            index <- which( mediatorNames[1,] == pathOrder[ 1, the_order_colNums1toN[j] ])

                          if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ the_order_colNums1toN[j] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)
                          }

                          if(mediator_model_type[index]=='lm'){

                                new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                                current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects
                          }

                    }


                    if(j == 1){
                                SAF[j] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                                previous_SAF <- SAF[j]
                    }else{
                                SAF[j] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)
                                # NB: The line "previous_SAF <- SAF[j]" has to come before the line "SAF[j] <- SAF[j] - previous_SAF"
                                previous_SAF_store <- SAF[j]

                                SAF[j] <- SAF[j] - previous_SAF

                                # NB: Have to define previous_SAF_store first and previous_SAF_store has to come before this line "SAF[j] <- SAF[j] - previous_SAF"
                                previous_SAF <- previous_SAF_store
                    }


                  }

###################################################################################################
###################################################################################################
# END OF NEW CODE TO ADD IN
###################################################################################################
###################################################################################################

    # MOC Don't need this line
    #SAF <- SAF/sum(w*sim_disease_current_population)
    SAF_mat[i,] <- SAF[reverse_order]
    order_mat[i,] <- the_order
    reverse_order_mat[i,] <- reverse_order
    if(i %% 100 == 0){
      flush.console()
      print(i)
    }
}

###################################################################################################
###################################################################################################
### CODE ADD IN FOR SEQUENTIAL AND JOINT PSPAF
###################################################################################################
###################################################################################################
###############################
###############################


    if(exact){
                # calculations are for joint PSPAFs rather than sequential PSPAFs
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
                  if( col_list[1:N][perm_mat[i,number_rf_new]] == riskfactor_col){
                         # assumes riskfactor is binary?
                         current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                         joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                  }else{
                            index <- which( mediatorNames[1,] == pathOrder[ 1, perm_mat[i,number_rf_new] ])

                            if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                # Set the risk factor to refval i.e. 0
                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                # # ## Add in to keep class of column after it is simulated
                                # keepColumnClass <- class(current_mat[,mediator_col[index]])
                                #
                                # temp_current_mat <- current_mat

                                # ## if is.factor(current_mat[,mediator_col[index]]) then do_sim outputs it as a character vector and impact_fraction outputs an error because the two classes for data and new_data differ as factor and character so need to keep class as character
                                # ## whereas it works if the class is numeric
                                # if( is.factor(current_mat[,mediator_col[index]]) ){
                                #
                                #   current_mat[,mediator_col[index]] <- factor(
                                #                                              do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,number_rf_new] ],
                                #                                              current_mat = current_mat_riskfactor_refval,
                                #                                              model = mediator_models[[index]],
                                #                                              SN=TRUE),
                                #                                              levels=levels(current_mat[,mediator_col[index]] )
                                #                                             )
                                #
                                # }else{
                                # # SET SN=TRUE as dont want SN=FALSE option
                                      current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,number_rf_new] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)
                                # }

                                # if( class(current_mat[,mediator_col[index]]) != keepColumnClass ){
                                #
                                #      if( keepColumnClass == ""factor)
                                #
                                #        for(i in 1:S){
                                #         if(is.factor(data[,i]) && length(levels(data[,i]))==2) data[,i] <- factor(as.numeric(data[,i]==levels(data[,i])[2]),levels=c(0,1))
                                #            if(!is.factor(data[,i]) && length(unique(data[,i]))==2 && is.numeric(data[,i])) data[,i] <- factor(as.numeric(data[,i]==max(data[,i])), levels=c(0,1))
                                #         if(is.character(data[,i])) data[,i] <- factor(data[,i],levels=sort(unique(data[,i])))
                                #
                                #       }
                                #
                                #       ## Add in to keep class of column after it is simulated
                                #       # class(current_mat[,mediator_col[index]]) <- class(temp_current_mat[,mediator_col[index]])
                                #
                                # }

                                joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)
                            }

                          if(mediator_model_type[index]=='lm'){

                              new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                              mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                              current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects

                              joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)
                          }

                  }

            }

            if(start_again==TRUE){

                  current_mat <- data
                  for(j in 1:number_rf_new){

                    if(col_list[1:N][perm_mat[i,j]] == riskfactor_col ){

                          current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)



                    }else{

                          index <- which( mediatorNames[1,] == pathOrder[ 1, perm_mat[i,j] ])

                          if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                ## if is.factor(current_mat[,mediator_col[index]]) then do_sim outputs it as a character vector and impact_fraction outputs an error because the two classes for data and new_data differ as factor and character so need to keep class as character
                                ## whereas it works if the class is numeric
                                # if( is.factor(current_mat[,mediator_col[index]]) ){
                                #
                                #       current_mat[,mediator_col[index]] <- factor(
                                #                                             do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,j] ],
                                #                                             current_mat = current_mat_riskfactor_refval,
                                #                                             model = mediator_models[[index]],
                                #                                             SN=TRUE),
                                #                                             levels=levels(current_mat[,mediator_col[index]] )
                                #                                             )
                                # }else{

                                      # MOC fix error here pathspecific_col_list_orig[ perm_mat[i,number_rf_new] ] should be pathspecific_col_list_orig[ perm_mat[i,j] ]
                                      current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,j] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)
                                # }

                          }

                          if(mediator_model_type[index]=='lm'){

                                new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                                current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects

                          }

                    }

                  }

                  joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

            }

    }
###################################################################################################
###################################################################################################
# END OF NEW CODE TO ADD IN
###################################################################################################
###################################################################################################


  }

  if(exact){
    # joint_PAF_vec <- joint_PAF_vec/sum(w*no_intervention)
    SAF_mat_exact <- matrix(0,nrow=N,ncol=N)
    rownames(SAF_mat_exact) <- paste('riskfactor ',1:N)
    colnames(SAF_mat_exact) <- paste('position ',1:N)
    for(i in 1:N){ # risk factor i
      for(j in 1:N){ # position j

        if(j < N) rows_to_look_at <- (1:nsim)[apply(perm_mat[,1:j,drop=FALSE],1,function(x){any(x==i)}) & perm_mat[,j]>0 & perm_mat[,j+1]==0]
        if(j == N) rows_to_look_at <- (1:nsim)[perm_mat[,N]>0]
        for(k in 1:length(rows_to_look_at)){
          joint_PAF_match_row <- 0
          if(j > 1){
          match_row <- perm_mat[rows_to_look_at[k],]
          match_row <- setdiff(match_row,i)
          match_row <- match_row[1:(j-1)]
          match_row <- (1:nsim)[apply(perm_mat,1,function(x){all(x[1:(j-1)]==match_row)&all(x[j:N]==0)})]
          joint_PAF_match_row <- joint_PAF_vec[match_row]
          }
          SAF_mat_exact[i,j] <- ((k-1)/k)*SAF_mat_exact[i,j]+(joint_PAF_vec[rows_to_look_at[k]]-joint_PAF_match_row)/k
          }
      }
    }

    #Each row has K! = Sum_{r=1 to K}((K-1) choose (r-1)) x (r-1)! x (K-r)!)  permutations.
      weights_exact <- numeric(N) # only used when exact calculation used
    for(i in 1:N){
      ## deleting the lines
      #   weights_exact[i] <- 1
      # if(i>1 && (N-i)>=1)
        weights_exact[i] <- choose(N-1, i-1)*factorial(i-1)*factorial(N-i)
    }

    average_PAF <- apply(SAF_mat_exact,1,function(x){weighted.mean(x, w=weights_exact)})
    SAF_mat_exact <- t(SAF_mat_exact)
    colnames(SAF_mat_exact) <- colnames(data)[col_list][1:N]
    names(average_PAF) <- colnames(data)[col_list][1:N]
    return(list(SAF_mat=SAF_mat_exact,average_PAF=average_PAF,joint_PAF=joint_PAF_vec[N]))
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

#' Calculation of average and sequential paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of models corresponding for the outcome variables in node_vec, with parents as described in parent_vec.  This list must be in the same order as node_vec and parent_list
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param exact logical.  Default TRUE. If TRUE, an efficient calculation is used to calculate average PAF, which enables the average PAF from N! permutations, over all N risk factors to be calculated with only 2^N-1 operations.  If FALSE, permutations are sampled
#' @param nsim  Default NULL Number of random permutations used to calculate average and sequential PAF.  If correct_order is set to an integer value, nsim is reset to the largest integer multiple of correct_order that is less than the number of permutations implied by correct_order.
#' @param correct_order Default 3.  This enforces stratified sampling of permutations where the first correct_order positions of the sampled permutations are evenly distributed over the integers 1 ... n, n being the number of risk factors of interest, over the sampled permutations.  The other positions are randomly sampled.  This automatically sets the number of simulations.  For interest, if n=10 and correct_order=3, nsim is set to factorial(n)/factorial(n-correct_order).  This special resampling reduces Monte Carlo variation in estimated average and sequential PAFs.
#' @params vars A subset of risk factors for which we want to calculate average, sequential and joint PAF
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendent nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE)
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available methods
#' @param ci_level Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the confidence interval (when ci=TRUE)
#' @param ci_level_ME Numeric.  Default 0.95. A number between 0 and 1 specifying the level of the margin of error for the point estimate
#' @return A dataframe with average joint and sequential PAF for all risk factors in node_vec (or alternatively a subset of those risk factors if specified in vars).
#'
#' @examples
#' set.seed(15122021)
#' library(dplyr)
#' library(devtools)
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=parallel::detectCores())
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure","early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,parent_alcohol,parent_stress,parent_high_blood_pressure,parent_lipids,parent_waist_hip_ratio,parent_early_stage_heart_disease,parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure","lipids","waist_hip_ratio","early_stage_heart_disease","diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list, node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)", spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' response_model <- glm(case ~ region*ns(age,df=5)+sex*ns(age,df=5) +  education + exercise + ns(waist_hip_ratio,df=5)+ smoking + alcohol + stress + high_blood_pressure + ns(lipids, knots = quantile(lipids,c(.25,0.5,0.75),na.rm=TRUE), Boundary.knots = quantile(lipids,c(.001,0.90),na.rm=TRUE))+ns(waist_hip_ratio,df=5)+diet,weights=weights, data=stroke_reduced,family='binomial')
#' mediator_models <- list(
#'   glm(formula = high_blood_pressure ~ region * ns(age, df = 5) + sex*ns(age, df = 5) + education + exercise + diet + smoking + alcohol + stress, family = "binomial", data = stroke_reduced, weights = weights),
#'   lm(formula = lipids ~ region * ns(age, df = 5) + sex * ns(age, df = 5) + education + exercise + diet + smoking + alcohol + stress, data = stroke_reduced, weights = weights),
#'   lm(formula = waist_hip_ratio ~ region * ns(age, df = 5) + sex * ns(age,df = 5) + education + exercise + diet + smoking + alcohol + stress,weights = weights,data=stroke_reduced))
#' ## No confidence interval, no bootstrap
#' # average_pspaf_no_CI(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = NULL, alpha=0.05,
#'                   # vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
#'                   # riskfactor = "exercise", refval=0, calculation_method = "D")
#' ## Bootstrap
#' average_pspaf(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = 3,
#'                  vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
#'                  riskfactor = "exercise", refval=0, calculation_method = "D", ci=TRUE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95)
#'
#' # joint_pspaf(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"),response_model = response_model, mediator_models = mediator_models, riskfactor = "exercise", refval=0, calculation_method = "D", ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95)
average_pspaf <- function(data, model_list, parent_list, node_vec, prev=.09, exact=TRUE, nsim=NULL, correct_order=2, vars=NULL,response_model, mediator_models, riskfactor, refval, calculation_method, ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95){
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

  if(!ci) return(average_pspaf_no_CI(data=data, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nsim=nsim, correct_order=correct_order, alpha=1-ci_level_ME,vars=vars,exact=exact, response_model=response_model, mediator_models=mediator_models, riskfactor=riskfactor, refval=refval, calculation_method = calculation_method))
  ## MOC: NEED TO UPDATE Parameters
  # CHECK IF NEED MORE PARAMETERS
  res <- boot::boot(data=data,statistic=average_pspaf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, nsim=nsim, correct_order=correct_order, vars=vars, exact=exact, response_model = response_model, mediator_models = mediator_models, riskfactor = riskfactor, refval=refval, calculation_method = calculation_method)
  if(is.null(vars)) vars <- node_vec[1:(length(node_vec)-1)]

      return(extract_ci(res=res,model_type='glm',t_vector=c(paste0(rep(node_vec[node_vec %in% vars],times=rep(length(vars),length(vars))),'_',rep(1:length(vars),length(vars))),paste0("Average PAF ", node_vec[node_vec %in% vars]),'JointPAF'),ci_level=ci_level,ci_type=ci_type,continuous=TRUE))

}

# CHECK PARAMETERS ADDED IN
average_pspaf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09, nsim=100, correct_order=3, vars=NULL, exact=TRUE, response_model = response_model, mediator_models = mediator_models,
                  riskfactor, refval, calculation_method, ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95){


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
      # return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
            if( is.factor(current_mat[,colnum ]) ){
                   return( factor( apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}),
                                   levels=levels(current_mat[,colnum ] ) ) )
            }else{
                   return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
                 }
    }
    # glm
    if(length(grep("glm",model$call))>0){

      probs <- predict(model,newdata=current_mat,type="response")
      if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
      # return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
            if( is.factor(current_mat[,colnum ]) ){
                     return( factor( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}),
                                     levels=levels(current_mat[,colnum ] ) ) )
              }else{
                     return( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}) )
                   }
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


  predict_df_discrete <- function(riskfactor, refval, data){

    if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

      stop("Riskfactor not in dataset.  Check spelling")

    }

    which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

    N <- nrow(data)
    riskfactor_vals <- data[,which_col]

    if(is.numeric(riskfactor_vals)){
      if(is.na(refval)) refval <- 0
      if(!all(riskfactor_vals %in% c(0,1)) || refval !=0){
        stop("Numeric risk factors must be 0/1, with the reference set to 0")

      }

    }

    if(is.numeric((data[,which_col]))) data[,which_col] <- rep(refval,N)
    if(is.character(data[,which_col])) data[,which_col] <- as.character(rep(refval,N))
    if(is.factor(data[,which_col])) data[,which_col] <- factor(rep(refval,N),levels = levels(data[,which_col]))

    return(data)

  }

  ##################################
  ########################  load in impact fraction functions:
  impact_fraction <- function(model, data, new_data, calculation_method="D",prev=NULL,ci=FALSE,boot_rep=100,t_vector=NULL, ci_level=0.95, ci_type=c("norm")){


    if(!is.data.frame(data)){
      stop(
        "data must be a dataframe object")
    }

    # remove data not used to fit model
    data <- data[row.names(data) %in% row.names(model.frame(model)),]
    new_data <- new_data[row.names(data) %in% row.names(model.frame(model)),]

    if(!is.data.frame(new_data)){
      stop(
        "new_data must be a dataframe object")
    }

    if(ncol(data)!=ncol(new_data) || nrow(data)!=nrow(new_data)){
      stop(
        "new_data must be the same dimensions as data")
    }

    if(!all(as.character(lapply(data,class))==as.character(lapply(new_data,class))) || !all(colnames(data)==colnames(new_data))){
      stop(
        "Data types and column names in new_data must match data types and column names in data.  To do this, try creating new_data from data")
    }

    if(!calculation_method %in% c("B","D")){
      stop(
        "Calculation of PAF only possible using the (B)ruzzi or (D)irect method.  Please supply either B or D")
    }
    response <- as.character(model$formula)[2]

    model_type <- NULL
    if(grepl("^glm$",as.character(model$call)[1],perl=TRUE)){
      model_type <- "glm"
      if (!as.character(model$family[1])=="binomial" & ! as.character(model$family[2]) %in% c("logit","log")) {
        stop(
          "The family must be binomial and link must be either log or logistic"
        )
      }
    }
    if(grepl("^coxph$",as.character(model$call)[1],perl=TRUE)){
      if("userCall" %in% names(model)){
        model_type <- "clogit"
        vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
        vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^strata\\((.*)\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^Surv\\(rep\\([0-9]*,[0-9]*\\),(.*)\\)$',replacement='\\1',x=vars)
        response <- vars[1]
      }else{
        model_type <- "coxph"
      }
    }
    if (is.null(model_type)) {
      stop(
        "Model must be either a glm, conditional logistic regression or Cox-proportional hazards regression"
      )
    }

    if (model_type=="coxph" && calculation_method=="B") {
      stop(
        "Bruzzi method unavailable with proportional hazards regression due to censoring.  Set method to direct instead"
      )
    }

    if(!is.null(prev) && model_type=="coxph"){

      stop(
        "Prevalence weighted estimation not appropriate for survival data sets"
      )


    }

    if(!is.double(t_vector) && !is.integer(t_vector) && model_type=="coxph"){

      stop(
        "Specify a numeric vector of times at which to calculate PAF"
      )


    }


    if(!is.null(prev) && (prev>=1 || prev<=0)){

      stop(
        "Prevalence should be a proportion (a number between 0 and 1)"
      )
    }

    N <- nrow(data)
    if(calculation_method=="B"){

      if(!ci) return(if_bruzzi(data, ind=1:N, model=model,model_type=model_type,new_data=new_data,response=response))
      if(ci){

        res <- boot::boot(data=data,statistic=if_bruzzi,R=boot_rep, model=model,model_type=model_type,new_data=new_data,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }
    if(calculation_method=="D"){

      if(!ci) return(if_direct(data,ind=1:N,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response))
      if(ci){
        res <- boot::boot(data=data,statistic=if_direct,R=boot_rep,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }

  }


  if_bruzzi <- function(data,ind, model,model_type,  new_data,response){

    N <- nrow(data)

    library(splines)
    library(survival)

    if(model_type == "clogit"){


      if(!all(ind==(1:N))){
        model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
        model_text <- paste0(model_text[2],model_text[1],model_text[3])
        strataname <- gsub(".*strata\\((.*)\\).*",replacement="\\1",x=model_text,perl=TRUE)

        # find strata variable
        strataids <- data[,colnames(data)==strataname]
        validids <- names(table(strataids))[table(strataids)==2]
        possibleresamples <- (1:nrow(data))[strataids %in% validids]
        ## order possible resamples indexes according to valid ids
        possibleresamples <- possibleresamples[order(strataids[strataids %in% validids])]
        totake <- sample(1:(length(possibleresamples)/2),length(possibleresamples)/2,replace=TRUE)
        resamples <- c(possibleresamples[2*totake],possibleresamples[2*totake-1])
        data <- data[resamples,]
        # avoid duplication of strata names
        data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))


        #refit model
        model_text <- paste0("clogit(",model_text,",data=data)")
        model <- eval(parse(text=model_text))

      }

      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #eta1 <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #eta2 <- the.mat%*%model$coefficients
      #oldRR <- exp(eta1)
      #newRR <- exp(eta2)
      oldRR <- predict(model,type="risk")
      newRR <- predict(model,type="risk",newdata=new_data)
      y <- data[,colnames(data)==response]
    }


    if(model_type == "glm"){

      if(!all(ind==(1:N))){

        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
        model <- eval(parse(text=model_text))

      }

      # predict on linear predictor scale
      oldRR <- exp(predict(model,newdata=data))
      newRR <- exp(predict(model,newdata=new_data))
      y <- data[,colnames(data)==response]
    }

    return(1 - mean(newRR[y==1]/oldRR[y==1]))

  }

  if_direct <- function(data, ind, model,model_type, new_data, prev,t_vector,response){

    library(splines)
    library(survival)

    N <- nrow(data)
    if(model_type == "coxph"){

      if(!all(ind==(1:N))){

        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("coxph(",model_text[2],",data=data)")
        model <- eval(parse(text=model_text))

      }
      cum_haz <- survival::basehaz(model, centered=FALSE)
      t_indices <- integer(length(t_vector))
      for(i in 1:length(t_vector)){
        t_indices[i] <- which.min(sapply(cum_haz[,2],function(x){(x-t_vector[i])^2}))
      }
      cum_haz <- cum_haz[t_indices,]
      oldhr <- predict(model,type="risk")
      newhr <- predict(model,newdata=new_data,type="risk")

      mean_probs_old <- 1 - apply(exp(-outer(cum_haz[,1],oldhr)),1,mean)

      mean_probs_new <- 1 - apply(exp(-outer(cum_haz[,1],newhr)),1,mean)

      PAF_vec <- (mean_probs_old - mean_probs_new)/mean_probs_old
      names(PAF_vec) <- paste0("t=",round(cum_haz[,2],2))
      return(PAF_vec)

    }



    add_term <- 0

    if(model_type=="clogit"){

      if(!all(ind==(1:N))){
        model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
        model_text <- paste0(model_text[2],model_text[1],model_text[3])
        strataname <- gsub(".*strata\\((.*)\\).*",replacement="\\1",x=model_text,perl=TRUE)

        # find strata variable
        strataids <- data[,colnames(data)==strataname]
        validids <- names(table(strataids))[table(strataids)==2]
        possibleresamples <- (1:nrow(data))[strataids %in% validids]
        ## order possible resamples indexes according to valid ids
        possibleresamples <- possibleresamples[order(strataids[strataids %in% validids])]
        totake <- sample(1:(length(possibleresamples)/2),length(possibleresamples)/2,replace=TRUE)
        resamples <- c(possibleresamples[2*totake],possibleresamples[2*totake-1])
        data <- data[resamples,]
        new_data <- new_data[resamples,]
        data <- data[resamples,]
        # avoid duplication of strata names
        data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))

        #refit model
        model_text <- paste0("clogit(",model_text,",data=data)")
        model <- eval(parse(text=model_text))

      }
      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #lp_old <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #lp_new <- the.mat%*%model$coefficients
      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model, newdata=new_data)
      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(!is.null(prev)){

        temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
        add_term <- uniroot(temp_fun, interval=c(-100,100))$root

      }
      probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
      probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))


    }else{
      # model is a glm

      if(!all(ind==(1:N))){
        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
        model <- eval(parse(text=model_text_u))

      }

      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model,newdata=new_data)

      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(as.character(model$family[2])=="logit"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
        probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))

      }
      if(as.character(model$family[2])=="log"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)
        probs_new <- exp(lp_new+add_term)

      }
    }
    return((sum(weights*probs_old)-sum(weights*probs_new))/sum(weights*probs_old))

  }
########


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
     nsim <- nrow(perm_mat)
       theorder <- apply(perm_mat,1,order_fun)
    perm_mat <- perm_mat[order(theorder,decreasing=FALSE),]
    }


  SAF_mat <- matrix(0,nrow=nsim,ncol=N)
  SAF_mat_2 <- matrix(0,nrow=nsim,ncol=N)
  order_mat <- matrix(0,nrow=nsim,ncol=N)
  reverse_order_mat <- matrix(0,nrow=nsim,ncol=N)
  joint_PAF_vec <- numeric(nsim) # only used when exact

  ######################################################################################
######################################################################################
######################################################################################
  ###################################################################################################
  ######
  ## Begininng of material to MOVE above for loop
  ######
  Num_forboot <- nrow(data)
  #Num <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- NULL
  if(grepl("^glm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "glm"
  if(grepl("^lm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "lm"
  if(grepl("^.*polr$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "polr"

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    if(grepl("^glm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "glm"
    if(grepl("^lm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "lm"
    if(grepl("^.*polr$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "polr"
  }

  ##### if(!all(ind==(1:N))){
  if(!all(ind==(1:Num_forboot))){

    data <- data[ind, ]
    if(response_model_type== "glm"){
      #browser()
      model_text <- as.character(response_model$call)
      if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"))")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"),weights=",model_text[5],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "lm"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "polr"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }
    for(i in 1:M){

      if(mediator_model_type[i]== "glm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"),weights=",model_text[5],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "lm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "polr"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }
    }



  }


mediatorNames <- as.data.frame( lapply(mediator_models,function(x){as.character(formula(x))[2]}))

pathOrder <- cbind(riskfactor,mediatorNames, as.character( formula(response_model) )[2])

pathspecific_col_list <- numeric(length(pathOrder))

N <- length(pathspecific_col_list)-1

for(i in 1:(N+1)) pathspecific_col_list[i] <- (1:ncol(data))[colnames(data)==pathOrder[1,i]]

  pathspecific_col_list_orig <- pathspecific_col_list

  riskfactor_col <- (1:length(colnames(data)))[colnames(data) %in% riskfactor]
  ######
  ## End of material to MOVE above for loop
  ######
######################################################################################
######################################################################################
######################################################################################

   for(i in 1:nsim){

  if(!exact){
    # if(is.null(correct_order)) the_order <- col_list[1:N][sample(1:N,N)]
    #### MOC CHECK
    if(is.null(correct_order)){
      the_order <- pathspecific_col_list[1:N][sample(1:N,N)]
      the_order_colNums1toN <- numeric(N)
      for(j in 1:N) the_order_colNums1toN[j] <- (1:N)[pathspecific_col_list_orig==the_order[j] ]
    }
    if(!is.null(correct_order)){

      the_order <- numeric(N)
      the_order[1:correct_order] <- perm_mat[i,1:correct_order]
      other_indexes <- setdiff(c(1:N),perm_mat[i,1:correct_order])
      ##### REPLACE WITH THIS  as.numeric( sample(as.character(other_indexes),N-correct_order) )
      # if(correct_order < N) the_order[(correct_order+1):N] <- sample(other_indexes,N-correct_order)
      if(correct_order < N) the_order[(correct_order+1):N] <- as.numeric( sample(as.character(other_indexes),N-correct_order) )
      the_order_colNums1toN <- the_order
      the_order <- col_list[1:N][the_order]
    }

    reverse_order <- numeric(N)
    for(j in 1:N) reverse_order[j] <- (1:N)[the_order==pathspecific_col_list[j]]

    current_mat <- data
    current_mat_2 <- data
    SAF <- numeric(N)
    SAF_2 <- numeric(N)
    no_intervention <- sim_disease_current_population

    # ############################################################
    # ############################################################
    # ### Previous Sequential PAF code
    # for(j in 1:length(the_order)){
    #
    #   current_mat <- sim_outnode(data,the_order[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
    #   current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
    #   SAF[j] <- (sum(w*no_intervention) - sum(w*current_mat[,col_list[N+1]]))
    #   no_intervention <- current_mat[,col_list[N+1]]
    #
    # }
    # ############################################################
    # ############################################################

###################################################################################################
###################################################################################################
### CODE ADD IN FOR Permutations  SEQUENTIAL AND JOINT PSPAF
###################################################################################################
###################################################################################################
###############################
###############################

                  # initialise previous_SAF at zero and same size as output from impact_fraction()
                  previous_SAF <- impact_fraction(model=response_model, data=data, new_data=data,calculation_method=calculation_method, prev=prev,ci=FALSE)
                  #current_mat <- data
                  # for(j in 1:number_rf_new){
                  for(j in 1:length(the_order)){

                    if( the_order[j]  == riskfactor_col ){

                          current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                    }else{
                            index <- which( mediatorNames[1,] == pathOrder[ 1, the_order_colNums1toN[j] ])

                          if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ the_order_colNums1toN[j] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)

                          }

                          if(mediator_model_type[index]=='lm'){

                                new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                                current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects
                          }

                    }


                    if(j == 1){
                                SAF[j] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                                previous_SAF <- SAF[j]
                    }else{
                                SAF[j] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                                # NB: The line "previous_SAF <- SAF[j]" has to come before the line "SAF[j] <- SAF[j] - previous_SAF"
                                previous_SAF_store <- SAF[j]

                                SAF[j] <- SAF[j] - previous_SAF

                                # NB: Have to define previous_SAF_store first and previous_SAF_store has to come before this line "SAF[j] <- SAF[j] - previous_SAF"
                                previous_SAF <- previous_SAF_store
                    }


                }

###################################################################################################
###################################################################################################
# END OF NEW CODE TO ADD IN
###################################################################################################
###################################################################################################

    # MOC Don't need this line
    # SAF <- SAF/sum(w*sim_disease_current_population)
    SAF_mat[i,] <- SAF[reverse_order]
    order_mat[i,] <- the_order
    reverse_order_mat[i,] <- reverse_order
  }


###################################################################################################
###################################################################################################
### CODE ADD IN FOR SEQUENTIAL AND JOINT PSPAF
###################################################################################################
###################################################################################################
###############################
###############################
    if(exact){
                # calculations are for joint PSPAFs rather than sequential PSPAFs
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
                  if( col_list[1:N][perm_mat[i,number_rf_new]] == riskfactor_col){
                         # assumes riskfactor is binary?
                         current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                         joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                  }else{
                            index <- which( mediatorNames[1,] == pathOrder[ 1, perm_mat[i,number_rf_new] ])

                            if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                # SET SN=TRUE as dont want SN=FALSE option
                                current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,number_rf_new] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)

                                joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)
                            }

                          if(mediator_model_type[index]=='lm'){

                              new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                              mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                              current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects

                              joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)


                          }

                  }

            }

            if(start_again==TRUE){

                  current_mat <- data
                  for(j in 1:number_rf_new){

                    if(col_list[1:N][perm_mat[i,j]] == riskfactor_col ){

                          current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                    }else{

                          index <- which( mediatorNames[1,] == pathOrder[ 1, perm_mat[i,j] ])

                          if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                current_mat[,mediator_col[index]] <- do_sim(colnum = pathspecific_col_list_orig[ perm_mat[i,j] ],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)
                          }

                          if(mediator_model_type[index]=='lm'){

                                new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                                current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects

                          }

                    }

                  }

                  joint_PAF_vec[i] <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

            }

    }
###################################################################################################
###################################################################################################
# END OF NEW CODE TO ADD IN
###################################################################################################
###################################################################################################

  }


  if(exact){
    # joint_PAF_vec <- joint_PAF_vec/sum(w*no_intervention)
    SAF_mat_exact <- matrix(0,nrow=N,ncol=N)
    rownames(SAF_mat_exact) <- paste('riskfactor ',1:N)
    colnames(SAF_mat_exact) <- paste('position ',1:N)
    for(i in 1:N){ # risk factor i
      for(j in 1:N){ # position j

        if(j < N) rows_to_look_at <- (1:nsim)[apply(perm_mat[,1:j,drop=FALSE],1,function(x){any(x==i)}) & perm_mat[,j]>0 & perm_mat[,j+1]==0]
        if(j == N) rows_to_look_at <- (1:nsim)[perm_mat[,N]>0]
        for(k in 1:length(rows_to_look_at)){
          joint_PAF_match_row <- 0
          if(j > 1){
            match_row <- perm_mat[rows_to_look_at[k],]
            match_row <- setdiff(match_row,i)
            match_row <- match_row[1:(j-1)]
            match_row <- (1:nsim)[apply(perm_mat,1,function(x){all(x[1:(j-1)]==match_row)&all(x[j:N]==0)})]
            joint_PAF_match_row <- joint_PAF_vec[match_row]
          }
          SAF_mat_exact[i,j] <- ((k-1)/k)*SAF_mat_exact[i,j]+(joint_PAF_vec[rows_to_look_at[k]]-joint_PAF_match_row)/k
        }
      }
    }

    # Each row has K! = Sum_{r=1 to K}((K-1) choose (r-1)) x (r-1)! x (K-r)!)  permutations.
    weights_exact <- numeric(N) # only used when exact calculation used
    for(i in 1:N){
      # weights_exact[i] <- 1
      # if(i>1 && (N-i)>=1)
        weights_exact[i] <- choose(N-1, i-1)*factorial(i-1)*factorial(N-i)
    }

    average_PAF <- apply(SAF_mat_exact,1,function(x){weighted.mean(x, w=weights_exact)})
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
    # return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
      if( is.factor(current_mat[,colnum ]) ){
                   return( factor( apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}),
                                                        levels=levels(current_mat[,colnum ] ) ) )
        }else{
                   return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
             }
  }
  # glm
  if(length(grep("glm",model$call))>0){

    probs <- predict(model,newdata=current_mat,type="response")
    if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
    # return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
    if( is.factor(current_mat[,colnum ]) ){
             return( factor( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}),
                                                                                      levels=levels(current_mat[,colnum ] ) ) )
      }else{
             return( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}) )
           }
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
#' Automatic fitting models for Bayesian network.
#'
#' Main effects models are fit by default.  For continuosu variables, lm is used, for binary (numeric 0/1 variables), glm is used and for factor valued variables polr is used.  For factors, ensure that the factor levels are ordered by increasing levels of risk.  If interactions are required for certain models, it is advisable to populate the elements of model_list separately.
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @param prev  Prevalence of disease.  Set to NULL for cohort or cross sectional studies
#' @param common character text for part of the model formula that doesn't involve any variable in node_vec.  Useful for specifying confounders involved in all models automatically
#' @param spline_nodes  Vector of continuous variable names that are fit as splines (when involved as parents)
#' @param df_spline_nodes How many df for each spline node
#' @return A list of models corresponding to node_vec and parent_vec.
#' @export
#'
#' @examples
#' More complicated example (slower to run)
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure","early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,parent_alcohol,parent_stress,parent_high_blood_pressure,parent_lipids,parent_waist_hip_ratio,parent_early_stage_heart_disease,parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure","lipids","waist_hip_ratio","early_stage_heart_disease","diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list, node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)", spline_nodes = c("waist_hip_ratio","lipids","diet"))
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

order_fun <- function(x){

  N <- length(x)
  sum <- 0
  for(i in 1:N){
    sum <- sum + x[i]*(N+1)^(N-i)
  }
  return(sum)
}


# Add in this function
predict_df_discrete <- function(riskfactor, refval, data){

    if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

      stop("Riskfactor not in dataset.  Check spelling")

    }

    which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

    N <- nrow(data)
    riskfactor_vals <- data[,which_col]

    if(is.numeric(riskfactor_vals)){
      if(is.na(refval)) refval <- 0
      if(!all(riskfactor_vals %in% c(0,1)) || refval !=0){
        stop("Numeric risk factors must be 0/1, with the reference set to 0")

      }

    }

    if(is.numeric((data[,which_col]))) data[,which_col] <- rep(refval,N)
    if(is.character(data[,which_col])) data[,which_col] <- as.character(rep(refval,N))
    if(is.factor(data[,which_col])) data[,which_col] <- factor(rep(refval,N),levels = levels(data[,which_col]))

    return(data)

}


##################################  the same functions as above are replicated here - but only return joint_PSPAF

#' Calculation of joint paf taking into account risk factor sequencing
#'
#' @param data Data frame. A dataframe containing variables used for fitting the models.  Must contain all variables used in fitting
#' @param model_list List.  A list of models corresponding for the outcome variables in node_vec, with parents as described in parent_vec.  This list must be in the same order as node_vec and parent_list
#' @param parent_list A list.  The ith element is the vector of variable names that are direct causes of ith variable in node_vec
#' @param node_vec A vector corresponding to the nodes in the Bayesian network.  This must be specified from root to leaves - that is ancestors in the causal graph for a particular node are positioned before their descendants.  If this condition is false the function will return an error.
#' @params vars A subset of risk factors for which we want to calculate average, sequential and joint PAF
#' @param ci Logical. If TRUE, a bootstrap confidence interval is computed along with a point estimate (default FALSE).  If ci=FALSE, only a point estimate is produced.  A simulation procedure (sampling permutations and also simulating the effects of eliminating risk factors over the descendent nodes in a Bayesian network) is required to produce the point estimates.  The point estimate will change on repated runs of the function.  The margin of error of the point estimate is given when ci=FALSE
#' @param boot_rep Integer.  Number of bootstrap replications (Only necessary to specify if ci=TRUE)
#' @param ci_type Character.  Default norm.  A vector specifying the types of confidence interval desired.  "norm", "basic", "perc" and "bca" are the available metho
#' @param ci_level Numeric.  Confidence level.  Default 0.95
#' @export
#'
#' @examples
#' set.seed(15122021)
#' library(dplyr)
#' library(devtools)
#' library(splines)
#' library(survival)
#' library(parallel)
#' options(boot.parallel="snow")
#' options(boot.ncpus=parallel::detectCores())
#' parent_exercise <- c("education")
#' parent_diet <- c("education")
#' parent_smoking <- c("education")
#' parent_alcohol <- c("education")
#' parent_stress <- c("education")
#' parent_high_blood_pressure <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_lipids <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_waist_hip_ratio <- c("education","exercise","diet","smoking","alcohol","stress")
#' parent_early_stage_heart_disease <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_diabetes <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure")
#' parent_case <- c("education","exercise","diet","smoking","alcohol","stress","lipids","waist_hip_ratio","high_blood_pressure","early_stage_heart_disease","diabetes")
#' parent_list <- list(parent_exercise,parent_diet,parent_smoking,parent_alcohol,parent_stress,parent_high_blood_pressure,parent_lipids,parent_waist_hip_ratio,parent_early_stage_heart_disease,parent_diabetes,parent_case)
#' node_vec=c("exercise","diet","smoking","alcohol","stress","high_blood_pressure","lipids","waist_hip_ratio","early_stage_heart_disease","diabetes","case")
#' model_list=automatic_fit(data=stroke_reduced, parent_list=parent_list, node_vec=node_vec, prev=.0035,common="region*ns(age,df=5)+sex*ns(age,df=5)", spline_nodes = c("waist_hip_ratio","lipids","diet"))
#' response_model <- glm(case ~ region*ns(age,df=5)+sex*ns(age,df=5) +  education + exercise + ns(waist_hip_ratio,df=5)+ smoking + alcohol + stress + high_blood_pressure + ns(lipids, knots = quantile(lipids,c(.25,0.5,0.75),na.rm=TRUE), Boundary.knots = quantile(lipids,c(.001,0.90),na.rm=TRUE))+ns(waist_hip_ratio,df=5)+diet,weights=weights, data=stroke_reduced,family='binomial')
#' mediator_models <- list(
#'   glm(formula = high_blood_pressure ~ region * ns(age, df = 5) + sex*ns(age, df = 5) + education + exercise + diet + smoking + alcohol + stress, family = "binomial", data = stroke_reduced, weights = weights),
#'   lm(formula = lipids ~ region * ns(age, df = 5) + sex * ns(age, df = 5) + education + exercise + diet + smoking + alcohol + stress, data = stroke_reduced, weights = weights),
#'   lm(formula = waist_hip_ratio ~ region * ns(age, df = 5) + sex * ns(age,df = 5) + education + exercise + diet + smoking + alcohol + stress,weights = weights,data=stroke_reduced))
#' ## No confidence interval, no bootstrap
#' # average_pspaf_no_CI(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = NULL, alpha=0.05,
#'                   # vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
#'                   # riskfactor = "exercise", refval=0, calculation_method = "D")
#' ## Bootstrap
#' # # average_pspaf(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = 3,
#'                 #  vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
#'                  # riskfactor = "exercise", refval=0, calculation_method = "D", ci=TRUE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95)
#'
#' joint_pspaf(data=stroke_reduced, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, vars = c("exercise","high_blood_pressure","lipids","waist_hip_ratio"),response_model = response_model, mediator_models = mediator_models, riskfactor = "exercise", refval=0, calculation_method = "D", ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95)
joint_pspaf <- function(data, model_list, parent_list, node_vec, prev=.09, vars=NULL,response_model, mediator_models, riskfactor, refval, calculation_method, ci=FALSE,boot_rep=100, ci_type=c("norm"),ci_level=0.95){
  if(!node_order(parent_list=parent_list,node_vec=node_vec)){
    stop("ancestors must be specified before descendants in node_vec")
  }
  if(!is.null(vars) & !all(vars %in% node_vec)){
    stop("Not all requested variables are in node_vec.  Check spelling")
  }
if(!ci) return(joint_pspaf_inner(data=data,ind=1:nrow(data), model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev,vars=vars, response_model=response_model, mediator_models=mediator_models, riskfactor=riskfactor, refval=refval, calculation_method=calculation_method))
  res <- boot::boot(data=data,statistic=joint_pspaf_inner,R=boot_rep,model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev=prev, vars=vars, response_model=response_model, mediator_models=mediator_models, riskfactor=riskfactor, refval=refval, calculation_method=calculation_method)
  return(boot::boot.ci(res,type=ci_type))

}


joint_pspaf_inner <- function(data, ind, model_list, parent_list, node_vec, prev=.09,vars=NULL, response_model, mediator_models, riskfactor, refval, calculation_method ){

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
      # return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
      if( is.factor(current_mat[,colnum ]) ){
                return( factor( apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}),
                                   levels=levels(current_mat[,colnum ] ) ) )
        }else{
                return(apply(probs,1,function(x){base::sample(mynames,size=1,prob=x)}))
             }

    }
    # glm
    if(length(grep("glm",model$call))>0){

      probs <- predict(model,newdata=current_mat,type="response")
      if(is.null(levels(current_mat[,colnum]))) return(apply(cbind(1-probs,probs),1,function(x){base::sample(c(0,1),size=1,prob=x)}))
      # return(apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}))
      if( is.factor(current_mat[,colnum ]) ){
                 return( factor( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}),
                                     levels=levels(current_mat[,colnum ] ) ) )
        }else{
                return( apply(cbind(1-probs,probs),1,function(x){base::sample(levels(current_mat[,colnum]),size=1,prob=x)}) )
             }
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
    ########################  load in impact fraction functions:
  impact_fraction <- function(model, data, new_data, calculation_method="B",prev=NULL,ci=FALSE,boot_rep=100,t_vector=NULL, ci_level=0.95, ci_type=c("norm")){


    if(!is.data.frame(data)){
      stop(
        "data must be a dataframe object")
    }

    # remove data not used to fit model
    data <- data[row.names(data) %in% row.names(model.frame(model)),]
    new_data <- new_data[row.names(data) %in% row.names(model.frame(model)),]

    if(!is.data.frame(new_data)){
      stop(
        "new_data must be a dataframe object")
    }

    if(ncol(data)!=ncol(new_data) || nrow(data)!=nrow(new_data)){
      stop(
        "new_data must be the same dimensions as data")
    }

    if(!all(as.character(lapply(data,class))==as.character(lapply(new_data,class))) || !all(colnames(data)==colnames(new_data))){
      stop(
        "Data types and column names in new_data must match data types and column names in data.  To do this, try creating new_data from data")
    }

    if(!calculation_method %in% c("B","D")){
      stop(
        "Calculation of PAF only possible using the (B)ruzzi or (D)irect method.  Please supply either B or D")
    }
    response <- as.character(model$formula)[2]

    model_type <- NULL
    if(grepl("^glm$",as.character(model$call)[1],perl=TRUE)){
      model_type <- "glm"
      if (!as.character(model$family[1])=="binomial" & ! as.character(model$family[2]) %in% c("logit","log")) {
        stop(
          "The family must be binomial and link must be either log or logistic"
        )
      }
    }
    if(grepl("^coxph$",as.character(model$call)[1],perl=TRUE)){
      if("userCall" %in% names(model)){
        model_type <- "clogit"
        vars <- gsub(pattern=' ',replacement='',x=unlist(strsplit(as.character(model$call)[2],split="[~*+]")))
        vars <- gsub(pattern='^ns\\((.*),df=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^ns\\((.*),knots=.*\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^strata\\((.*)\\)$',replacement='\\1',x=vars)
        vars <- gsub(pattern='^Surv\\(rep\\([0-9]*,[0-9]*\\),(.*)\\)$',replacement='\\1',x=vars)
        response <- vars[1]
      }else{
        model_type <- "coxph"
      }
    }
    if (is.null(model_type)) {
      stop(
        "Model must be either a glm, conditional logistic regression or Cox-proportional hazards regression"
      )
    }

    if (model_type=="coxph" && calculation_method=="B") {
      stop(
        "Bruzzi method unavailable with proportional hazards regression due to censoring.  Set method to direct instead"
      )
    }

    if(!is.null(prev) && model_type=="coxph"){

      stop(
        "Prevalence weighted estimation not appropriate for survival data sets"
      )


    }

    if(!is.double(t_vector) && !is.integer(t_vector) && model_type=="coxph"){

      stop(
        "Specify a numeric vector of times at which to calculate PAF"
      )


    }


    if(!is.null(prev) && (prev>=1 || prev<=0)){

      stop(
        "Prevalence should be a proportion (a number between 0 and 1)"
      )
    }

    N <- nrow(data)
    if(calculation_method=="B"){

      if(!ci) return(if_bruzzi(data, ind=1:N, model=model,model_type=model_type,new_data=new_data,response=response))
      if(ci){

        res <- boot::boot(data=data,statistic=if_bruzzi,R=boot_rep, model=model,model_type=model_type,new_data=new_data,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }
    if(calculation_method=="D"){

      if(!ci) return(if_direct(data,ind=1:N,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response))
      if(ci){
        res <- boot::boot(data=data,statistic=if_direct,R=boot_rep,model=model, model_type=model_type, new_data=new_data, prev=prev,t_vector=t_vector,response=response)
        return(extract_ci(res=res,model_type=model_type,t_vector=t_vector,ci_level=ci_level,ci_type=ci_type))
      }
    }

  }


  if_bruzzi <- function(data,ind, model,model_type,  new_data,response){

    N <- nrow(data)

    library(splines)
    library(survival)

    if(model_type == "clogit"){


      if(!all(ind==(1:N))){
        model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
        model_text <- paste0(model_text[2],model_text[1],model_text[3])
        strataname <- gsub(".*strata\\((.*)\\).*",replacement="\\1",x=model_text,perl=TRUE)

        # find strata variable
        strataids <- data[,colnames(data)==strataname]
        validids <- names(table(strataids))[table(strataids)==2]
        possibleresamples <- (1:nrow(data))[strataids %in% validids]
        ## order possible resamples indexes according to valid ids
        possibleresamples <- possibleresamples[order(strataids[strataids %in% validids])]
        totake <- sample(1:(length(possibleresamples)/2),length(possibleresamples)/2,replace=TRUE)
        resamples <- c(possibleresamples[2*totake],possibleresamples[2*totake-1])
        data <- data[resamples,]
        # avoid duplication of strata names
        data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))


        #refit model
        model_text <- paste0("clogit(",model_text,",data=data)")
        model <- eval(parse(text=model_text))

      }

      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #eta1 <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #eta2 <- the.mat%*%model$coefficients
      #oldRR <- exp(eta1)
      #newRR <- exp(eta2)
      oldRR <- predict(model,type="risk")
      newRR <- predict(model,type="risk",newdata=new_data)
      y <- data[,colnames(data)==response]
    }


    if(model_type == "glm"){

      if(!all(ind==(1:N))){

        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
        model <- eval(parse(text=model_text))

      }

      # predict on linear predictor scale
      oldRR <- exp(predict(model,newdata=data))
      newRR <- exp(predict(model,newdata=new_data))
      y <- data[,colnames(data)==response]
    }

    return(1 - mean(newRR[y==1]/oldRR[y==1]))

  }

  if_direct <- function(data, ind, model,model_type, new_data, prev,t_vector,response){

    library(splines)
    library(survival)

    N <- nrow(data)
    if(model_type == "coxph"){

      if(!all(ind==(1:N))){

        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        model_text <- paste0("coxph(",model_text[2],",data=data)")
        model <- eval(parse(text=model_text))

      }
      cum_haz <- survival::basehaz(model, centered=FALSE)
      t_indices <- integer(length(t_vector))
      for(i in 1:length(t_vector)){
        t_indices[i] <- which.min(sapply(cum_haz[,2],function(x){(x-t_vector[i])^2}))
      }
      cum_haz <- cum_haz[t_indices,]
      oldhr <- predict(model,type="risk")
      newhr <- predict(model,newdata=new_data,type="risk")

      mean_probs_old <- 1 - apply(exp(-outer(cum_haz[,1],oldhr)),1,mean)

      mean_probs_new <- 1 - apply(exp(-outer(cum_haz[,1],newhr)),1,mean)

      PAF_vec <- (mean_probs_old - mean_probs_new)/mean_probs_old
      names(PAF_vec) <- paste0("t=",round(cum_haz[,2],2))
      return(PAF_vec)

    }



    add_term <- 0

    if(model_type=="clogit"){

      if(!all(ind==(1:N))){
        model_text <- as.character(eval(parse(text=as.character(model$userCall)[2])))
        model_text <- paste0(model_text[2],model_text[1],model_text[3])
        strataname <- gsub(".*strata\\((.*)\\).*",replacement="\\1",x=model_text,perl=TRUE)

        # find strata variable
        strataids <- data[,colnames(data)==strataname]
        validids <- names(table(strataids))[table(strataids)==2]
        possibleresamples <- (1:nrow(data))[strataids %in% validids]
        ## order possible resamples indexes according to valid ids
        possibleresamples <- possibleresamples[order(strataids[strataids %in% validids])]
        totake <- sample(1:(length(possibleresamples)/2),length(possibleresamples)/2,replace=TRUE)
        resamples <- c(possibleresamples[2*totake],possibleresamples[2*totake-1])
        data <- data[resamples,]
        new_data <- new_data[resamples,]
        data <- data[resamples,]
        # avoid duplication of strata names
        data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))
        new_data <- new_data[resamples,]
        # avoid duplication of strata names
        new_data[,colnames(data)==strataname] <- c(1:length(totake),1:length(totake))

        #refit model
        model_text <- paste0("clogit(",model_text,",data=data)")
        model <- eval(parse(text=model_text))

      }
      # clogit inherits predictions from coxph.  They seem strange at first but are equivalent to predictions from the following code which takes longer to run and so is commented out
      #model$coefficients[is.na(model$coefficients)] <- 0
      #them <- model.matrix(model)
      #them[is.na(model.matrix(model))] <- 0
      #lp_old <- them%*%model$coefficients
      #the.mat <- model.matrix(as.formula(paste("~",gsub(paste('+ strata(',strataname,')',sep=''),'',x=as.character(model$formula[3]),fixed=TRUE),sep="")),data=new_data)
      #the.mat <- the.mat[,-1,drop=FALSE]
      #the.mat[is.na(the.mat)] <- 0
      #lp_new <- the.mat%*%model$coefficients
      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model, newdata=new_data)
      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(!is.null(prev)){

        temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
        add_term <- uniroot(temp_fun, interval=c(-100,100))$root

      }
      probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
      probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))


    }else{
      # model is a glm

      if(!all(ind==(1:N))){
        data <- data[ind, ]
        new_data <- new_data[ind, ]
        model_text <- as.character(model$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(model)[2]),"),weights=",model_text[5],")")
        model <- eval(parse(text=model_text_u))

      }

      lp_old <- predict(model,newdata=data)
      lp_new <- predict(model,newdata=new_data)

      y <- data[,colnames(data)==response]
      N <- nrow(data)
      weights <- rep(1, N)
      if(!is.null(prev)){

        data_prev <- mean(as.numeric(y==1))
        weights[y==0] <- (1-prev)/(1-data_prev)
        weights[y==1] <- prev/data_prev

      }


      if(as.character(model$family[2])=="logit"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old)/(1+exp(c+lp_old)),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)/(1+exp(lp_old+add_term))
        probs_new <- exp(lp_new+add_term)/(1+exp(lp_new+add_term))

      }
      if(as.character(model$family[2])=="log"){

        if(!is.null(prev)){

          temp_fun <- function(c){weighted.mean(exp(c+lp_old),w=weights)-prev}
          add_term <- uniroot(temp_fun, interval=c(-100,100))$root

        }
        probs_old <- exp(lp_old+add_term)
        probs_new <- exp(lp_new+add_term)

      }
    }
    return((sum(weights*probs_old)-sum(weights*probs_new))/sum(weights*probs_old))

  }

  predict_df_discrete <- function(riskfactor, refval, data){

    if(all(!grepl(paste0("^",riskfactor,"$"),colnames(data),perl=TRUE))){

      stop("Riskfactor not in dataset.  Check spelling")

    }

    which_col <- grep(paste0("^",riskfactor,"$"),colnames(data))

    N <- nrow(data)
    riskfactor_vals <- data[,which_col]

    if(is.numeric(riskfactor_vals)){
      if(is.na(refval)) refval <- 0
      if(!all(riskfactor_vals %in% c(0,1)) || refval !=0){
        stop("Numeric risk factors must be 0/1, with the reference set to 0")

      }

    }

    if(is.numeric((data[,which_col]))) data[,which_col] <- rep(refval,N)
    if(is.character(data[,which_col])) data[,which_col] <- as.character(rep(refval,N))
    if(is.factor(data[,which_col])) data[,which_col] <- factor(rep(refval,N),levels = levels(data[,which_col]))

    return(data)

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

  ######################################################################################
######################################################################################
######################################################################################
  ###################################################################################################
  ######
  ## Begininng of material to MOVE above for loop
  ######
  Num_forboot <- nrow(data)
  #Num <- nrow(data)
  riskfactor_col <- grep(paste0('^',riskfactor,'$'),colnames(data),perl=TRUE)
  M <- length(mediator_models)
  mediator_col <- rep(1,M)
  for(i in 1:M) mediator_col[i] <- grep(as.character(formula(mediator_models[[i]]))[2],colnames(data),perl=TRUE)

  response_model_type <- NULL
  if(grepl("^glm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "glm"
  if(grepl("^lm$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "lm"
  if(grepl("^.*polr$",as.character(response_model$call)[1],perl=TRUE)) response_model_type <- "polr"

  mediator_model_type <- rep(NULL, M)
  for(i in 1:M){
    if(grepl("^glm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "glm"
    if(grepl("^lm$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "lm"
    if(grepl("^.*polr$",as.character(mediator_models[[i]]$call)[1],perl=TRUE)) mediator_model_type[i] <- "polr"
  }

  ##### if(!all(ind==(1:N))){
  if(!all(ind==(1:Num_forboot))){

    data <- data[ind, ]
    if(response_model_type== "glm"){
      #browser()
      model_text <- as.character(response_model$call)
      if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"))")
      if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(response_model)[2]),"),weights=",model_text[5],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "lm"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }

    if(response_model_type == "polr"){
      model_text <- as.character(response_model$call)
      if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
      if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
      response_model <- eval(parse(text=model_text_u))
    }
    for(i in 1:M){

      if(mediator_model_type[i]== "glm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==4) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"))")
        if(length(model_text)==5) model_text_u <- paste0("glm(",model_text[2],",data=data, family=binomial(link=",as.character(family(mediator_models[[i]])[2]),"),weights=",model_text[5],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "lm"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("lm(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("lm(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }

      if(mediator_model_type[i] == "polr"){
        model_text <- as.character(mediator_models[[i]]$call)
        if(length(model_text)==3) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data)")
        if(length(model_text)==4) model_text_u <- paste0("MASS::polr(",model_text[2],",data=data, weights=",model_text[4],")")
        mediator_models[[i]] <- eval(parse(text=model_text_u))
      }
    }



  }


mediatorNames <- as.data.frame( lapply(mediator_models,function(x){as.character(formula(x))[2]}))

pathOrder <- cbind(riskfactor,mediatorNames, as.character( formula(response_model) )[2])

pathspecific_col_list <- numeric(length(pathOrder))

N <- length(pathspecific_col_list)-1

for(i in 1:(N+1)) pathspecific_col_list[i] <- (1:ncol(data))[colnames(data)==pathOrder[1,i]]

  pathspecific_col_list_orig <- pathspecific_col_list

  riskfactor_col <- (1:length(colnames(data)))[colnames(data) %in% riskfactor]
  ######
  ## End of material to MOVE above for loop
  ######
######################################################################################
######################################################################################
######################################################################################


#### Previous Code
# current_mat <- data
#       for(j in 1:length(col_list)){
#
#         current_mat <- sim_outnode(data,col_list[j],current_mat,parent_list=parent_list,col_list=col_list_orig,model_list=model_list)
#         current_mat[,col_list[N+1]] <- predict(model_list[[length(node_vec)]],newdata=current_mat,type="response")
#
#       }
#       return(jointPAF=(sum(w*sim_disease_current_population)-sum(w*current_mat[,col_list[N+1]]))/sum(w*sim_disease_current_population))

  ######################################################################################
  current_mat <- data
                    for(j in 1:(length(col_list) - 1) ){

                    if(col_list[j] == riskfactor_col ){

                          current_mat <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                    }else{

                          index <- which( mediator_col == col_list[j]  )

                          if(mediator_model_type[index]=='glm' || mediator_model_type[index]=='polr'){

                                current_mat_riskfactor_refval <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                current_mat[,col_list[j]] <- do_sim(colnum = col_list[j],
                                                                            current_mat = current_mat_riskfactor_refval,
                                                                            model = mediator_models[[index]],
                                                                            SN=TRUE)
                          }

                          if(mediator_model_type[index]=='lm'){

                                new_data_direct <- predict_df_discrete(riskfactor=riskfactor, refval=refval, data = current_mat)

                                mediator_effects <- predict(mediator_models[[index]]) - predict(mediator_models[[index]],new_data_direct)

                                current_mat[,mediator_col[index]] <- current_mat[,mediator_col[index]] - mediator_effects

                          }

                    }

                  }

                  joint_PAF_vec <- impact_fraction(model=response_model, data=data, new_data=current_mat,calculation_method=calculation_method, prev=prev,ci=FALSE)

                  return(jointPAF= joint_PAF_vec )
########################################

}
