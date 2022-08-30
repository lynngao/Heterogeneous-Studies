library(caret)
library(pROC)
library(operators)
library(readr)
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", "nnls",  "glmnet", 
             "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", "parallel", "doParallel",  
             "ranger", "scales"), require, character.only=TRUE))
#change to your directory for helper file
source("/Users/lynngao/Desktop/simulation of heterogenity/helper.R")


command_args <- commandArgs(trailingOnly=TRUE)
#command_args <- c(1,1,100,10,10,1.075,"logit",
                  #"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_yu_ctr.rds", 
                  #"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_feng_ctr.rds")
if(length(command_args)!=10){stop("Not enough input parameters!")} 
#population difference alpha(0,0.5,1), library size factor lambda(0.5,1,2), sample size, 
#genes, disease effect factor(1,2,5), dataset1, dataset2, saved table location
alpha = as.numeric(command_args[1]) #population difference factor
lambda = as.numeric(command_args[2]) #library size factor
sample_size = as.numeric(command_args[3])
num_gene = as.numeric(command_args[4])
overlap_gene = as.numeric(command_args[5])
ed = as.numeric(command_args[6])
method = as.character(command_args[7])

count_data1 <- readRDS(command_args[8]) #Asian population
count_data2 <- readRDS(command_args[9]) #Caucacian population

#get fixed library size as median of all samples
library_size = 1000000

top_otu1 = colnames(count_data1)[order(apply(count_data1, 2, var), decreasing = TRUE)[1:1000]]
count_data1 = count_data1[,colnames(count_data1)%in%top_otu1]
count_data1 = count_data1[,sort(colnames(count_data1))]
for (i in 1:length(top_otu1)){
  if (top_otu1[i] %!in% colnames(count_data2)){
    count_data2 = cbind(count_data2, 0)
    colnames(count_data2)[ncol(count_data2)] = top_otu1[i]  
  }
}
count_data2 = count_data2[,colnames(count_data2)%in%top_otu1]
count_data2 = count_data2[,sort(colnames(count_data2))]

#real probabilties
prob_data1 = as.vector(colMeans(count_data1)/sum(colMeans(count_data1)))
prob_data2 = as.vector(colMeans(count_data2)/sum(colMeans(count_data2)))

otu_list1 = which(prob_data1!=0)
otu_list2 = which(prob_data2!=0)
set.seed(1)
otu_selected1 = sample(otu_list1, num_gene)
set.seed(1)
otu_selected2 = sample(otu_selected1, overlap_gene)


#pseudo probabitlies for controls
pseudo_prob <- function(alpha, prob1, prob2, num_gene, ed, otu_selected1, otu_selected2){
  otu_enrich = otu_selected1[1:(num_gene/2)]
  otu_deplete = otu_selected1[(num_gene/2+1):num_gene]
  
  prob_v1_ctr = alpha*prob1 + (1-alpha)*prob2
  prob_v2_ctr = prob2
  
  prob_v1_case = prob_v1_ctr
  prob_v1_case[otu_enrich] = prob_v1_ctr[otu_enrich]*ed
  prob_v1_case[otu_deplete] = prob_v1_ctr[otu_deplete]/ed
  
  prob_v2_case = prob_v2_ctr
  for (i in 1:length(otu_selected2)){
    if (otu_selected2[i] %in% otu_enrich){
      prob_v2_case[otu_selected2[i]] = prob_v2_ctr[otu_selected2[i]]*ed
    }
    if (otu_selected2[i] %in% otu_deplete){
      prob_v2_case[otu_selected2[i]] = prob_v2_ctr[otu_selected2[i]]/ed
    }
  }
  
  prob_v1_case = prob_v1_case/sum(prob_v1_case)
  prob_v2_case = prob_v2_case/sum(prob_v2_case)
  return(list(prob_v1_ctr, prob_v2_ctr, prob_v1_case, prob_v2_case))
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

prob_v = pseudo_prob(alpha, prob_data1, prob_data2, num_gene, ed, otu_selected1, otu_selected2)
prob_v1_ctr = prob_v[[1]]
prob_v2_ctr = prob_v[[2]]
prob_v1_case = prob_v[[3]]
prob_v2_case = prob_v[[4]]


# Set a dataframe to save results
perf_df <- as.data.frame(matrix(NA, nrow=102, ncol=30))
colnames(perf_df) = c("Training1","Training2", "Training1_combat","Training2_combat","Merged", "ComBat_Merged",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_auc", "LOSO_auc", "rank_mean","rank_geo_mean","rank_Stuart", "rank_RRA","rank_Bayesian", 
                      "Avg_combat", "n_Avg_combat", "CS_Avg_combat", "Reg_a_combat", "Reg_s_combat", "val_auc_combat", "LOSO_auc_combat", "rank_mean_combat","rank_geo_mean_combat","rank_Stuart_combat", "rank_RRA_combat","rank_Bayesian_combat")
rownames(perf_df) = c(1:100, "mean", "sd")

for (i in 1:100) {
  simulation1 = as.data.frame(generate_samples(alpha, lambda, prob_v1_ctr, prob_v1_case, count_data1, sample_size)) #training set1
  simulation2 = as.data.frame(generate_samples(alpha, lambda, prob_v1_ctr, prob_v1_case, count_data1, sample_size)) #training set2
  simulation3 = as.data.frame(generate_samples(alpha, lambda, prob_v2_ctr, prob_v2_case, count_data2, sample_size)) #testing set
  
  ## ComBat normalize two training sets
  simulation1_combat = combat_norm(simulation3, simulation1)
  simulation2_combat = combat_norm(simulation3, simulation2)
  
  ## Merge two training sets into one
  merged_training = rbind(simulation1, simulation2)
  merged_training_combat = rbind(simulation1_combat, simulation2_combat)
  
  ## Split testing set into validation and testing
  testing = simulation3
  sample_case = sample.int(n=nrow(testing[testing$status == "case",]),size=floor(0.5*nrow(testing[testing$status == "case",])),replace=F)
  sample_ctr = sample.int(n=nrow(testing[testing$status == "ctr",]),size=floor(0.5*nrow(testing[testing$status == "ctr",])),replace=F)
  val_case = testing[testing$status == "case",][sample_case,]
  val_ctr = testing[testing$status == "ctr",][sample_ctr,]
  val = rbind(val_case,val_ctr)
  testing = testing[rownames(testing)%!in%rownames(val),]
  
  ## Obtain predictions from learner trained within each training set
  ml1 = ml_model(simulation1, method)
  ml2 = ml_model(simulation2, method)
  ml1_combat = ml_model(simulation1_combat, method)
  ml2_combat = ml_model(simulation2_combat, method)
  pred_prob1 = predict(ml1, testing, type="prob")$case
  pred_prob2 = predict(ml2, testing, type="prob")$case
  pred_sgbatch_res <- list(pred_prob1, pred_prob2)
  pred_prob1_combat = predict(ml1_combat, testing, type="prob")$case
  pred_prob2_combat = predict(ml2_combat, testing, type="prob")$case
  pred_sgbatch_res_combat  <- list(pred_prob1_combat, pred_prob2_combat )
  
  ## Prediction from combine training together (Merged)
  ml_merged = ml_model(merged_training, method)
  pred_merged_res = predict(ml_merged, testing, type="prob")$case
  
  ## Prediction from training after batch adjustment (Merged combat)
  ml_merged_combat = ml_model(merged_training_combat, method)
  pred_combat_res = predict(ml_merged_combat, testing, type="prob")$case
  
  ## Prediction on validation dataset
  pred_sgbatch_res_val <- list(predict(ml1, val, type="prob")$case, predict(ml2, val, type="prob")$case)
  names(pred_sgbatch_res_val) <- paste0("Simulation", 1:2)
  pred_sgbatch_res_val_combat <- list(predict(ml1_combat, val, type="prob")$case, predict(ml2_combat, val, type="prob")$case)
  names(pred_sgbatch_res_val_combat) <- paste0("Simulation", 1:2)
  
  ##  Aggregate with different weights
  pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp)})
  pred_mat <- do.call(cbind, pred_test_lst)
  pred_test_lst_combat <- lapply(pred_sgbatch_res_combat, function(tmp){return(tmp)})
  pred_mat_combat <- do.call(cbind, pred_test_lst_combat)
  
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  pred_avg_combat <- rowMeans(pred_mat_combat)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% (as.matrix(rbind(nrow(simulation1), nrow(simulation2))) / (nrow(simulation1) + nrow(simulation2)))
  pred_N_avg_combat <- pred_mat_combat %*% (as.matrix(rbind(nrow(simulation1_combat), nrow(simulation2_combat))) / (nrow(simulation1_combat) + nrow(simulation2_combat)))
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=2, training = list(simulation1, simulation2), perf_name="mxe", method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_combat <- CS_zmatrix(n_batch=2, training = list(simulation1_combat, simulation2_combat), perf_name="mxe", method)
  cs_weights_seq_combat <- CS_weight(cs_zmat_combat)
  pred_cs_avg_combat <- pred_mat_combat %*% cs_weights_seq_combat
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=2, training = list(simulation1, simulation2), method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(list(simulation1, simulation2), nrow))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_combat <- Reg_SSL_pred(n_batch=2, training = list(simulation1_combat, simulation2_combat), method)
  reg_a_beta_combat <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_combat$coef), 
                                    n_seq=sapply(list(simulation1_combat, simulation2_combat), nrow))
  pred_reg_a_combat <- pred_mat_combat %*% reg_a_beta_combat
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = c(simulation1$status, simulation2$status)
  stacked_label = ifelse(stacked_label == "ctr", 0, 1)
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_combat <- do.call(rbind, reg_ssl_res_combat$pred)
  reg_s_beta_combat <- nnls(A=stacked_pred_combat, b=as.numeric(stacked_label))$x
  reg_s_beta_combat <- reg_s_beta_combat / sum(reg_s_beta_combat)
  pred_reg_s_combat <- pred_mat_combat %*% reg_s_beta_combat
  
  
  # Weighting by validation set performance
  coef1 = auc(val$status, pred_sgbatch_res_val$Simulation1)[1]
  coef2 = auc(val$status, pred_sgbatch_res_val$Simulation2)[1]
  pred_val_auc = as.matrix(coef1*pred_mat[,1] + coef2*pred_mat[,2])/(coef1+coef2)
  
  coef1_combat = auc(val$status, pred_sgbatch_res_val_combat$Simulation1)[1]
  coef2_combat = auc(val$status, pred_sgbatch_res_val_combat$Simulation2)[1]
  pred_val_auc_combat = as.matrix(coef1_combat*pred_mat_combat[,1] + coef2_combat*pred_mat_combat[,2])/(coef1_combat+coef2_combat)
  
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = 4, nrow = nrow(testing)))
  colnames(LOSO) = c("auc_ds1","auc_ds2","auc_ds1_combat","auc_ds2_combat")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = auc(testing$status[-k], pred_prob1[-k])
    LOSO[k,2] = auc(testing$status[-k], pred_prob2[-k])
    LOSO[k,3] = auc(testing$status[-k], pred_prob1_combat[-k])
    LOSO[k,4] = auc(testing$status[-k], pred_prob2_combat[-k])
  }
  LOSO[LOSO < 0.5] = 0.5
  pred_LOSO = ifelse(LOSO[,1]+LOSO[,2] == 1, 0.5, ((LOSO[,1]-0.5)*pred_prob1 + (LOSO[,2]-0.5)*pred_prob2)/ (LOSO[,1]+LOSO[,2]-1))
  pred_LOSO_combat = ifelse(LOSO[,3]+LOSO[,4] == 1, 0.5, ((LOSO[,3]-0.5)*pred_prob1_combat + (LOSO[,4]-0.5)*pred_prob2_combat)/ (LOSO[,3]+LOSO[,4]-1))
  
  
  # rank methods
  rank_mat = rank_method(testing, pred_prob1, pred_prob2)
  rank_mat_combat = rank_method(testing, pred_prob1_combat, pred_prob2_combat)
  
  
  ####  Evaluate performance
  tst_scores <- c(list(Training1 = pred_prob1, Training2 = pred_prob2, Training1_combat = pred_prob1_combat, Training2_combat = pred_prob2_combat, 
                       Merged=pred_merged_res, ComBat_Merged=pred_combat_res, 
                       Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO, 
                       rank_mean = rank_mat$mean, rank_geo_mean = rank_mat$geo_mean, rank_Stuart=rank_mat$Q, rank_RRA=rank_mat$rho, rank_Bayesian=rank_mat$Bayesian,
                       Avg_combat=pred_avg_combat, n_Avg_combat=pred_N_avg_combat, CS_Avg_combat=pred_cs_avg_combat, Reg_a_combat=pred_reg_a_combat, Reg_s_combat=pred_reg_s_combat, val_auc_combat = pred_val_auc_combat, LOSO_auc_combat = pred_LOSO_combat, 
                       rank_mean_combat = rank_mat_combat$mean, rank_geo_mean_combat = rank_mat_combat$geo_mean, rank_Stuart_combat=rank_mat_combat$Q, rank_RRA_combat=rank_mat_combat$rho, rank_Bayesian_combat=rank_mat_combat$Bayesian))
  
  perf_df[i,] <- sapply(tst_scores, function(preds){
    testing_auc <- auc(testing$status, as.vector(preds))
    return(testing_auc[1])})
}

# Get mean and sd
perf_df[101,] = colMeans(perf_df[1:100,])
perf_df[102,] = apply(perf_df[1:100,], 2, sd)

####  Output results
file_name <- sprintf('scenario3_%s__auc_a%s_l%s_s%s_g%s_e%s_overlap%s_test%s_logRelAbun.csv', 
                     gsub('.', '', method, fixed=T),
                     gsub('.', '', alpha, fixed=T), gsub('.', '', lambda, fixed=T), 
                     gsub('.', '', sample_size, fixed=T), gsub('.', '', num_gene, fixed=T), gsub('.', '', ed, fixed=T), gsub('.', '', overlap_gene, fixed=T),
                     gsub('.', '', sub(".*_count_*(.*?)_ctr.*", "\\1", command_args[9]), fixed=T))
write.csv(perf_df, paste0(command_args[10],file_name))
