library(caret)
library(pROC)
library(operators)
library(readr)
rm(list=ls())
setwd("/Users/lynngao/Desktop/simulation of heterogenity/scenario2_results/")
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", "nnls",  "glmnet", 
             "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", "parallel", "doParallel",  
             "ranger", "scales"), require, character.only=TRUE))
#change to your directory for helper file
source("/Users/lynngao/Desktop/simulation of heterogenity/helper.R")

####  Simulate data from same background population
command_args <- commandArgs(trailingOnly=TRUE)
#command_args <- c(0,1,300,100,10,1.075,"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_yu_ctr.rds",
                  #2,100,100,3,1,"logit")
alpha = as.numeric(command_args[1]) #population size difference
lambda = as.numeric(command_args[2]) #library size factor
sample_size1 = as.numeric(command_args[3])
sample_size2 = as.numeric(command_args[4])
num_gene = as.numeric(command_args[5])
ed = as.numeric(command_args[6])
count_data1 <- readRDS(command_args[7]) #Asian population


top_otu1 = colnames(count_data1)[order(apply(count_data1, 2, var), decreasing = TRUE)[1:1000]]
count_data1 = count_data1[,colnames(count_data1)%in%top_otu1]
count_data1 = count_data1[,sort(colnames(count_data1))]

#get fixed library size as median of all samples
library_size = 1000000

#real probabilties
prob_data1 = as.vector(colMeans(count_data1)/sum(colMeans(count_data1)))

otu_list = which(prob_data1!=0)
set.seed(1)
otu_selected = sample(otu_list, num_gene)

#pseudo probabitlies
pseudo_prob <- function(alpha, prob1, num_gene, ed, otu_selected){
  otu_enrich = otu_selected[1:(num_gene/2)]
  otu_deplete = otu_selected[(num_gene/2+1):num_gene]
  
  prob_v1_ctr = prob1
  prob_v1_case = prob_v1_ctr
  prob_v1_case[otu_enrich] = prob_v1_ctr[otu_enrich]*ed
  prob_v1_case[otu_deplete] = prob_v1_ctr[otu_deplete]/ed
  prob_v1_case = prob_v1_case/sum(prob_v1_case)
   
  return(list(prob_v1_ctr, prob_v1_case))
}

generate_samples <- function(alpha, lambda, prob_ctr, prob_case, count_data, sample_size){
  simulation_ctr = rmultinom(sample_size*0.5, size = library_size*lambda, prob = prob_ctr)
  simulation_case = rmultinom(sample_size*0.5, size = library_size*lambda, prob = prob_case)
  simulation_adjust = as.data.frame(t(cbind(simulation_ctr, simulation_case)))
  colnames(simulation_adjust) = colnames(count_data)
  simulation_adjust = simulation_adjust/rowSums(simulation_adjust) #change to relative abundance
  min = min(apply(simulation_adjust[,1:ncol(simulation_adjust)], 1, function(x) min(x[x>0])))
  simulation_adjust[simulation_adjust == 0] = min*0.65
  simulation_adjust = log(simulation_adjust)
  simulation_adjust$status = c(rep("ctr",sample_size*0.5), rep("case",sample_size*0.5))
  return(simulation_adjust)
}

prob_v = pseudo_prob(alpha, prob_data1, num_gene, ed, otu_selected)
prob_v1_ctr = prob_v[[1]]
prob_v1_case = prob_v[[2]]


####  Parameters 
#command_args <- commandArgs(trailingOnly=TRUE)  
#command_args <- c("20", "5", "1")
#if(length(command_args)!=3){stop("Not enough input parameters!")}

## Degree of batch effect (strength of signal)
N_batch <- as.numeric(command_args[8])
N_sample_size <- c(as.numeric(command_args[9]),as.numeric(command_args[10]))
max_batch_mean <- as.numeric(command_args[11])
max_batch_var <- as.numeric(command_args[12])
#N_sample_size <- as.numeric(command_args[1])   # number of samples per batch (case + control)
#max_batch_mean <- as.numeric(command_args[2]) 
# median 0.12, range(0, 12), recommended values 0 0.3 0.5 (additive)
#max_batch_var <- as.numeric(command_args[3]) 
# median 0.02, range(0.002, 0.33), recommended values 1 2 4 (multiplicative)
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha, 
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  
cat("\nBatch changes\n");
print(hyper_pars$hyper_mu);
print(hyper_pars$hyper_beta/(hyper_pars$hyper_alpha-1))
#print((hyper_pars$hyper_beta^2)/((hyper_pars$hyper_alpha-1)^2*(hyper_pars$hyper_alpha-2)))

## Pipeline
iterations <- 100
#exp_name <- sprintf('batchN%s_m%s_v%s', N_sample_size, 
#gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T)) 
exp_name <- sprintf('batch_m%s_v%s', 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T))
perf_measures <- "auc"    
method = as.character(command_args[13])
l_type <- method

####  Run pipeline
#c; #l_type = learner_types[2]
#start_time <- Sys.time()
for(ID in 1:iterations){
  simulation1 = as.data.frame(generate_samples(alpha, lambda, prob_v1_ctr, prob_v1_case, count_data1, sample_size1)) #training set1
  simulation2 = as.data.frame(generate_samples(alpha, lambda, prob_v1_ctr, prob_v1_case, count_data1, sample_size2)) #testing set
  
  # divide test data into training, test and validation
  training = simulation1
  sample_case = subset(simulation2, status=="case")
  sample_ctr = subset(simulation2, status=="ctr")
  val_case = sort(as.integer(sample(rownames(sample_case), size=floor(0.5*nrow(sample_case)), replace =F)))
  val_ctr = sort(as.integer(sample(rownames(sample_ctr), size=floor(0.5*nrow(sample_ctr)), replace =F)))
  val_id = sort(c(val_case,val_ctr))
  val = simulation2[rownames(simulation2)%in%val_id,]
  testing = simulation2[rownames(simulation2)%!in%val_id,]
  
  tst_scores_modlst <- cs_zmat_lst <- list()  # list of results for each model
  
  ## Subset training set in batches
  #batches_ind <- subsetBatch(condition=y_train, N_sample_size=N_sample_size, N_batch=N_batch)
  sample_case1 = sort(as.integer(sample(rownames(training[training$status=="case",]), size=0.5*N_sample_size[1], replace =F)))
  sample_ctr1 = sort(as.integer(sample(rownames(training[training$status=="ctr",]), size=0.5*N_sample_size[1], replace =F)))
  sample_id1 = sort(c(sample_case1,sample_ctr1))
  sample_case2 = sort(as.integer(sample(rownames(training[training$status=="case" & rownames(training)%!in%sample_id1,]), size=0.5*N_sample_size[2], replace =F)))
  sample_ctr2 = sort(as.integer(sample(rownames(training[training$status=="ctr" & rownames(training)%!in%sample_id1,]), size=0.5*N_sample_size[2], replace =F)))
  sample_id2 = sort(c(sample_case2,sample_ctr2))

  batch1 = training[rownames(training)%in%sample_id1,]
  batch2 = training[rownames(training)%in%sample_id2,]

  batches_ind = list(sample_id1, sample_id2)
  batch = c(rep(1,N_sample_size[1]), rep(2,N_sample_size[2]))
  
  ## Remove genes with only 0 values in any batch in current training set
  g1 = colnames(batch1[colSums(batch1[,-ncol(batch1)]) != 0])
  g2 = colnames(batch2[colSums(batch2[,-ncol(batch2)]) != 0])
  g_keep <- intersect(g1, g2)
  curr_test_expr <- testing[,g_keep]
  curr_val_expr <- val[,g_keep]
  batch1 = batch1[,g_keep]
  batch2 = batch2[,g_keep]
  curr_train_expr = rbind(batch1, batch2)
  
  ## Simulate batch effect 
  sim_batch_res <- simBatch(curr_train_expr, N_sample_size, batches_ind, batch, hyper_pars)
  #if(ID==1){print(round(do.call(cbind, sim_batch_res$batch_par), 3))}
  train_expr_batch <- sim_batch_res$new_dat
  batch1_simulated = train_expr_batch[1:N_sample_size[1],]
  batch2_simulated = train_expr_batch[(N_sample_size[1]+1):nrow(train_expr_batch),]

  ##ComBat normalization on batch and testing dataset (simulation2 as reference)
  batch1_combat <- combat_norm(simulation2, batch1_simulated)
  rownames(batch1_combat) = rownames(batch1_simulated)
  batch2_combat <- combat_norm(simulation2, batch2_simulated)
  rownames(batch2_combat) = rownames(batch2_simulated)
  
  train_expr_combat = rbind(batch1_combat, batch2_combat)
  train_expr_combat_norm <- train_expr_combat
  
  train_expr_norm <- as.data.frame(curr_train_expr)
  test_expr_norm <- curr_test_expr
  val_expr_norm <- curr_val_expr
  train_expr_batch_whole_norm <- train_expr_batch_norm <- as.data.frame(train_expr_batch)
  
  train_lst <- lapply(batches_ind, function(ind){train_expr_batch_norm[as.character(ind),]})
  train_lst_combat <- lapply(batches_ind, function(ind){train_expr_combat_norm[as.character(ind),]})
  
  ####  Training with single model
  print(sprintf("Simulation: %s, Model: %s", ID, l_type))
  
  
  ## Prediction from original training to test, without batch effect
  ml_merged <- ml_model(train_expr_norm, method)
  pred_base_res <- predict(ml_merged, test_expr_norm, type="prob")$case
  
  
  ## Prediction from training WITH batch effect to test
  ml_merged_batch <- ml_model(train_expr_batch_whole_norm, method)
  pred_batch_res <- predict(ml_merged_batch, test_expr_norm, type="prob")$case
  
  
  ##  Prediction from training after batch adjustment (Merged)
  ml_merged_batch_combat <- ml_model(train_expr_combat_norm, method)
  pred_combat_res <- predict(ml_merged_batch_combat, test_expr_norm, type="prob")$case
  
  
  ## Obtain predictions from learner trained within each batch
  ml1 <- ml_model(batch1_simulated, method)
  ml2 <- ml_model(batch2_simulated, method)
  ml1_combat <- ml_model(batch1_combat, method)
  ml2_combat <- ml_model(batch2_combat, method)
  pred_prob1 = predict(ml1, test_expr_norm, type="prob")$case
  pred_prob2 = predict(ml2, test_expr_norm, type="prob")$case
  pred_prob1_combat = predict(ml1_combat, test_expr_norm, type="prob")$case
  pred_prob2_combat = predict(ml2_combat, test_expr_norm, type="prob")$case
  pred_sgbatch_res <- list(pred_prob1, pred_prob2)
  pred_sgbatch_res_combat <- list(pred_prob1_combat, pred_prob2_combat)
  names(pred_sgbatch_res) <- paste0("Batch", 1:N_batch)
  names(pred_sgbatch_res_combat) <- paste0("Batch", 1:N_batch, "_combat")
  
  
  ##prediction on validation dataset
  pred_prob1_val = predict(ml1, val_expr_norm, type="prob")$case
  pred_prob2_val = predict(ml2, val_expr_norm, type="prob")$case
  pred_prob1_val_combat = predict(ml1_combat, val_expr_norm, type="prob")$case
  pred_prob2_val_combat = predict(ml2_combat, val_expr_norm, type="prob")$case
  pred_sgbatch_res_val <- list(pred_prob1_val, pred_prob2_val)
  names(pred_sgbatch_res_val) <- paste0("Batch", 1:N_batch)
  pred_sgbatch_res_val_combat <- list(pred_prob1_val_combat, pred_prob2_val_combat)
  names(pred_sgbatch_res_val_combat) <- paste0("Batch", 1:N_batch, "_combat")
  
  
  ##  Aggregate with different weights
  pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp)})
  pred_mat <- do.call(cbind, pred_test_lst)
  pred_test_lst_combat <- lapply(pred_sgbatch_res_combat, function(tmp){return(tmp)})
  pred_mat_combat <- do.call(cbind, pred_test_lst_combat)
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  pred_avg_combat <- rowMeans(pred_mat_combat)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / nrow(curr_train_expr))
  pred_N_avg_combat <- pred_mat_combat %*% (as.matrix(sapply(batches_ind, length)) / nrow(curr_train_expr))
  
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=N_batch, training = train_lst, perf_name="mxe", method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_combat <- CS_zmatrix(n_batch=N_batch, training = train_lst_combat, perf_name="mxe", method)
  cs_weights_seq_combat <- CS_weight(cs_zmat_combat)
  pred_cs_avg_combat <- pred_mat_combat %*% cs_weights_seq_combat
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=N_batch, training = train_lst, method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(batches_ind, length))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_combat <- Reg_SSL_pred(n_batch=N_batch, training = train_lst_combat, method)
  reg_a_beta_combat <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_combat$coef), 
                             n_seq=sapply(batches_ind, length))
  pred_reg_a_combat <- pred_mat_combat %*% reg_a_beta_combat
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = curr_train_expr$status
  stacked_label = ifelse(stacked_label == "ctr", 0, 1)
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_combat <- do.call(rbind, reg_ssl_res_combat$pred)
  reg_s_beta_combat <- nnls(A=stacked_pred_combat, b=as.numeric(stacked_label))$x
  reg_s_beta_combat <- reg_s_beta_combat / sum(reg_s_beta_combat)
  pred_reg_s_combat <- pred_mat_combat %*% reg_s_beta_combat

  
  # Weighting by validation set performance
  coef1 = auc(val_expr_norm$status, pred_sgbatch_res_val$Batch1)[1]
  coef2 = auc(val_expr_norm$status, pred_sgbatch_res_val$Batch2)[1]
  pred_val_auc = as.matrix(coef1*pred_mat[,1] + coef2*pred_mat[,2])/(coef1+coef2)
  
  coef1_combat = auc(val_expr_norm$status, pred_sgbatch_res_val_combat$Batch1)[1]
  coef2_combat = auc(val_expr_norm$status, pred_sgbatch_res_val_combat$Batch2)[1]
  pred_val_auc_combat = as.matrix(coef1_combat*pred_mat_combat[,1] + coef2_combat*pred_mat_combat[,2])/(coef1_combat+coef2_combat)


  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = N_batch*2, nrow = nrow(test_expr_norm)))
  colnames(LOSO) = c("auc_ds1","auc_ds2", "auc_ds1_combat","auc_ds2_combat")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = auc(testing$status[-k], pred_sgbatch_res$Batch1[-k])[1]
    LOSO[k,2] = auc(testing$status[-k], pred_sgbatch_res$Batch2[-k])[1]
    LOSO[k,3] = auc(testing$status[-k], pred_sgbatch_res_combat$Batch1[-k])[1]
    LOSO[k,4] = auc(testing$status[-k], pred_sgbatch_res_combat$Batch2[-k])[1]
  }
  LOSO[LOSO < 0.5] = 0.5
  pred_LOSO = as.matrix(((LOSO[,1]-0.5)*pred_sgbatch_res$Batch1 + (LOSO[,2]-0.5)*pred_sgbatch_res$Batch2) / (LOSO[,1]+LOSO[,2]-1))
  pred_LOSO_combat = as.matrix(((LOSO[,3]-0.5)*pred_sgbatch_res_combat$Batch1 + (LOSO[,4]-0.5)*pred_sgbatch_res_combat$Batch2) / (LOSO[,3]+LOSO[,4]-1))
  for (j in 1:length(pred_LOSO)){
    if(is.nan(pred_LOSO[j])){
      pred_LOSO[j] = 0.5
    }
    if(is.nan(pred_LOSO_combat[j])){
      pred_LOSO_combat[j] = 0.5
    }
  }
  
  
  # rank methods
  rank_mat = rank_method(test_expr_norm, pred_prob1, pred_prob2)
  rank_mat_combat = rank_method(test_expr_norm, pred_prob1_combat, pred_prob2_combat)
  
  
  ####  Evaluate performance
  tst_scores <- c(list(NoBatch=pred_base_res, Batch=pred_batch_res), 
                  pred_test_lst, pred_test_lst_combat,
                  list(ComBat_Merged=pred_combat_res, Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO,
                       rank_mean = rank_mat$mean, rank_geo_mean = rank_mat$geo_mean, rank_Stuart=rank_mat$Q, rank_RRA=rank_mat$rho, rank_Bayesian=rank_mat$Bayesian, 
                       Avg_combat=pred_avg_combat, n_Avg_combat=pred_N_avg_combat, CS_Avg_combat=pred_cs_avg_combat, Reg_a_combat=pred_reg_a_combat, Reg_s_combat=pred_reg_s_combat, val_auc_combat = pred_val_auc_combat, LOSO_auc_combat = pred_LOSO_combat,
                       rank_mean_combat = rank_mat_combat$mean, rank_geo_mean_combat = rank_mat_combat$geo_mean, rank_Stuart_combat=rank_mat_combat$Q, rank_RRA_combat=rank_mat_combat$rho, rank_Bayesian_combat=rank_mat_combat$Bayesian))
  
  perf_df <- lapply(perf_measures, function(perf_name){
    as.data.frame(t(sapply(tst_scores, function(preds){
      testing_auc <- auc(test_expr_norm$status, as.vector(preds))
      return(testing_auc[1])
    })))
  })
  names(perf_df) <- perf_measures
  perf_df <- do.call(rbind, perf_df)
  
  ####  Output results
  first_file <- !file.exists(sprintf('/home1/yilingao/proj/simulation_heterogeneity/results/integration/scenario2_results/%s_auc_%s.csv', l_type, exp_name))
  wo <- sapply(perf_measures, function(perf_name){
    write.table(perf_df[perf_name, ], sprintf('/home1/yilingao/proj/simulation_heterogeneity/results/integration/scenario2_results/%s_%s_%s.csv', l_type, perf_name, exp_name),
                append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
  })
  tst_scores_modlst[[l_type]] <- tst_scores
  cs_zmat_lst[[l_type]] <- cs_zmat
}
#end_time <- Sys.time()
#print(end_time - start_time)

