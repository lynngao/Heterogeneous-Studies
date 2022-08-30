library(caret)
library(pROC)
library(operators)
library(readr)
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", "nnls",  "glmnet", 
             "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", "parallel", "doParallel",  
             "ranger", "scales"), require, character.only=TRUE))
source("/Users/lynngao/Desktop/simulation of heterogenity/helper.R")

command_args <- commandArgs(trailingOnly=TRUE)
#command_args <- c("/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/A.rds", 
                  #"/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/C.rds", 
                  #"/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/D.rds", 
                  #"/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/E.rds", 
                  #"/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/F.rds", 
                  #"/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/G.rds", "rf")
if(length(command_args)!=8){stop("Not enough input parameters!")} 


count_data1 <- readRDS(command_args[1]) #training 1
count_data2 <- readRDS(command_args[2]) #training 2
count_data3 <- readRDS(command_args[3]) #training 3
count_data4 <- readRDS(command_args[4]) #training 4
count_data5 <- readRDS(command_args[5]) #training 5
count_data6 <- readRDS(command_args[6]) #testing

method = as.character(command_args[8])

get_count_data <- function(count_data){
  status_data = count_data$status
  count_data <- count_data[,-ncol(count_data)]
  top_otu = colnames(count_data)[order(apply(count_data, 2, var), decreasing = TRUE)[1:1000]]
  count_data = count_data[,colnames(count_data)%in%top_otu]
  return(list(status_data, as.data.frame(count_data)))
}

temp1 = get_count_data(count_data1)
temp2 = get_count_data(count_data2)
temp3 = get_count_data(count_data3)
temp4 = get_count_data(count_data4)
temp5 = get_count_data(count_data5)
temp6 = get_count_data(count_data6)
status_data1 = temp1[[1]]
count_data1 = temp1[[2]]
status_data2 = temp2[[1]]
count_data2 = temp2[[2]]
status_data3 = temp3[[1]]
count_data3 = temp3[[2]]
status_data4 = temp4[[1]]
count_data4 = temp4[[2]]
status_data5 = temp5[[1]]
count_data5 = temp5[[2]]
status_data6 = temp6[[1]]
count_data6 = temp6[[2]]

#get union of species
union_species = unique(c(colnames(count_data1), colnames(count_data2),colnames(count_data3),colnames(count_data4),colnames(count_data5)))

trim_count_data <- function(union_species, count_data){
  for (i in 1:length(union_species)){
    if (union_species[i] %!in% colnames(count_data)){
      count_data = cbind(count_data, 0)
      colnames(count_data)[ncol(count_data)] = union_species[i]  
    }
  }
  count_data = count_data[,colnames(count_data)%in%union_species]
  count_data = count_data[,sort(colnames(count_data))]
  return(as.data.frame(count_data))
}

count_data1 = trim_count_data(union_species, count_data1)
count_data2 = trim_count_data(union_species, count_data2)
count_data3 = trim_count_data(union_species, count_data3)
count_data4 = trim_count_data(union_species, count_data4)
count_data5 = trim_count_data(union_species, count_data5)
count_data6 = trim_count_data(union_species, count_data6)

get_rel_abun <- function(ds, status){
  #ds = ds/rowSums(ds) #change to relative abundance
  #min = min(apply(ds[,1:ncol(ds)], 1, function(x) min(x[x>0])))
  #ds[ds == 0] = min*0.65
  #ds = log(ds)
  ds$status = status
  return(as.data.frame(ds))
}

ds1 = get_rel_abun(count_data1,status_data1)
ds2 = get_rel_abun(count_data2,status_data2)
ds3 = get_rel_abun(count_data3,status_data3)
ds4 = get_rel_abun(count_data4,status_data4)
ds5 = get_rel_abun(count_data5,status_data5)
ds6 = get_rel_abun(count_data6,status_data6)

# Set a dataframe to save results
perf_df <- as.data.frame(matrix(NA, nrow=32, ncol=36))
colnames(perf_df) = c("Training1","Training2","Training3","Training4","Training5",
                      "Training1_combat","Training2_combat","Training3_combat","Training4_combat","Training5_combat",
                      "Merged", "ComBat_Merged",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_auc", "LOSO_auc", "rank_mean","rank_geo_mean","rank_Stuart", "rank_RRA","rank_Bayesian", 
                      "Avg_combat", "n_Avg_combat", "CS_Avg_combat", "Reg_a_combat", "Reg_s_combat", "val_auc_combat", "LOSO_auc_combat", "rank_mean_combat","rank_geo_mean_combat","rank_Stuart_combat", "rank_RRA_combat","rank_Bayesian_combat")
rownames(perf_df) = c(1:30, "mean", "sd")

for (i in 1:30) {
  ## ComBat normalize training sets
  ds1_combat = combat_norm(ds6, ds1)
  ds2_combat = combat_norm(ds6, ds2)
  ds3_combat = combat_norm(ds6, ds3)
  ds4_combat = combat_norm(ds6, ds4)
  ds5_combat = combat_norm(ds6, ds5)
  
  ## Merge two training sets into one
  merged_training = rbind(ds1, ds2, ds3, ds4, ds5)
  merged_training_combat = rbind(ds1_combat, ds2_combat, ds3_combat, ds4_combat, ds5_combat)
  
  ## Split testing set into validation and testing
  testing = ds6
  sample_case = sample.int(n=nrow(testing[testing$status == "case",]),size=floor(0.5*nrow(testing[testing$status == "case",])),replace=F)
  sample_ctr = sample.int(n=nrow(testing[testing$status == "ctr",]),size=floor(0.5*nrow(testing[testing$status == "ctr",])),replace=F)
  val_case = testing[testing$status == "case",][sample_case,]
  val_ctr = testing[testing$status == "ctr",][sample_ctr,]
  val = rbind(val_case,val_ctr)
  testing = testing[rownames(testing)%!in%rownames(val),]
  
  ## Obtain predictions from learner trained within each training set
  rf1 = ml_model(ds1, method)
  rf2 = ml_model(ds2, method)
  rf3 = ml_model(ds3, method)
  rf4 = ml_model(ds4, method)
  rf5 = ml_model(ds5, method)
  rf1_combat = ml_model(ds1_combat, method)
  rf2_combat = ml_model(ds2_combat, method)
  rf3_combat = ml_model(ds3_combat, method)
  rf4_combat = ml_model(ds4_combat, method)
  rf5_combat = ml_model(ds5_combat, method)
  
  pred_prob1 = predict(rf1, testing, type="prob")$case
  pred_prob2 = predict(rf2, testing, type="prob")$case
  pred_prob3 = predict(rf3, testing, type="prob")$case
  pred_prob4 = predict(rf4, testing, type="prob")$case
  pred_prob5 = predict(rf5, testing, type="prob")$case
  pred_sgbatch_res <- list(pred_prob1, pred_prob2, pred_prob3, pred_prob4, pred_prob5)
  
  pred_prob1_combat = predict(rf1_combat, testing, type="prob")$case
  pred_prob2_combat = predict(rf2_combat, testing, type="prob")$case
  pred_prob3_combat = predict(rf3_combat, testing, type="prob")$case
  pred_prob4_combat = predict(rf4_combat, testing, type="prob")$case
  pred_prob5_combat = predict(rf5_combat, testing, type="prob")$case
  pred_sgbatch_res_combat  <- list(pred_prob1_combat, pred_prob2_combat, pred_prob3_combat, pred_prob4_combat, pred_prob5_combat)
  
  ## Prediction from combine training together (Merged)
  rf_merged = ml_model(merged_training, method)
  pred_merged_res = predict(rf_merged, testing, type="prob")$case
  
  ## Prediction from training after batch adjustment (Merged combat)
  rf_merged_combat = ml_model(merged_training_combat, method)
  pred_combat_res = predict(rf_merged_combat, testing, type="prob")$case
  
  ## Prediction on validation dataset
  pred_sgbatch_res_val <- list(predict(rf1, val, type="prob")$case, predict(rf2, val, type="prob")$case,
                               predict(rf3, val, type="prob")$case, predict(rf4, val, type="prob")$case,
                               predict(rf5, val, type="prob")$case)
  names(pred_sgbatch_res_val) <- paste0("ds", 1:5)
  pred_sgbatch_res_val_combat <- list(predict(rf1_combat, val, type="prob")$case, predict(rf2_combat, val, type="prob")$case,
                                      predict(rf3_combat, val, type="prob")$case, predict(rf4_combat, val, type="prob")$case,
                                      predict(rf5_combat, val, type="prob")$case)
  names(pred_sgbatch_res_val_combat) <- paste0("ds", 1:5)
  
  ##  Aggregate with different weights
  pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp)})
  pred_mat <- do.call(cbind, pred_test_lst)
  pred_test_lst_combat <- lapply(pred_sgbatch_res_combat, function(tmp){return(tmp)})
  pred_mat_combat <- do.call(cbind, pred_test_lst_combat)
  
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  pred_avg_combat <- rowMeans(pred_mat_combat)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% (as.matrix(rbind(nrow(ds1), nrow(ds2), nrow(ds3), nrow(ds4), nrow(ds5))) / (nrow(ds1) + nrow(ds2) + nrow(ds3) + nrow(ds4) + nrow(ds5)))
  pred_N_avg_combat <- pred_mat_combat %*% (as.matrix(rbind(nrow(ds1_combat), nrow(ds2_combat), nrow(ds3_combat), nrow(ds4_combat), nrow(ds5_combat))) / (nrow(ds1_combat) + nrow(ds2_combat) + nrow(ds3_combat) + nrow(ds4_combat) + nrow(ds5_combat)))
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=5, training = list(ds1, ds2, ds3, ds4, ds5), perf_name="mxe")
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_combat <- CS_zmatrix(n_batch=5, training = list(ds1_combat, ds2_combat, ds3_combat, ds4_combat, ds5_combat), perf_name="mxe")
  cs_weights_seq_combat <- CS_weight(cs_zmat_combat)
  pred_cs_avg_combat <- pred_mat_combat %*% cs_weights_seq_combat
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=5, training = list(ds1, ds2, ds3, ds4, ds5))
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(list(ds1, ds2, ds3, ds4, ds5), nrow))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_combat <- Reg_SSL_pred(n_batch=5, training = list(ds1_combat, ds2_combat, ds3_combat, ds4_combat, ds5_combat))
  reg_a_beta_combat <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_combat$coef), 
                                    n_seq=sapply(list(ds1_combat, ds2_combat, ds3_combat, ds4_combat, ds5_combat), nrow))
  pred_reg_a_combat <- pred_mat_combat %*% reg_a_beta_combat
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = c(ds1$status, ds2$status, ds3$status, ds4$status, ds5$status)
  stacked_label = ifelse(stacked_label == "ctr", 0, 1)
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_combat <- do.call(rbind, reg_ssl_res_combat$pred)
  reg_s_beta_combat <- nnls(A=stacked_pred_combat, b=as.numeric(stacked_label))$x
  reg_s_beta_combat <- reg_s_beta_combat / sum(reg_s_beta_combat)
  pred_reg_s_combat <- pred_mat_combat %*% reg_s_beta_combat
  
  
  # Weighting by validation set performance
  coef1 = auc(val$status, pred_sgbatch_res_val$ds1)[1]
  coef2 = auc(val$status, pred_sgbatch_res_val$ds2)[1]
  coef3 = auc(val$status, pred_sgbatch_res_val$ds3)[1]
  coef4 = auc(val$status, pred_sgbatch_res_val$ds4)[1]
  coef5 = auc(val$status, pred_sgbatch_res_val$ds5)[1]
  pred_val_auc = as.matrix(coef1*pred_mat[,1] + coef2*pred_mat[,2] + coef3*pred_mat[,3] + coef4*pred_mat[,4] + coef5*pred_mat[,5])/(coef1+coef2+coef3+coef4+coef5)
  
  coef1_combat = auc(val$status, pred_sgbatch_res_val_combat$ds1)[1]
  coef2_combat = auc(val$status, pred_sgbatch_res_val_combat$ds2)[1]
  coef3_combat = auc(val$status, pred_sgbatch_res_val_combat$ds3)[1]
  coef4_combat = auc(val$status, pred_sgbatch_res_val_combat$ds4)[1]
  coef5_combat = auc(val$status, pred_sgbatch_res_val_combat$ds5)[1]
  pred_val_auc_combat = as.matrix(coef1_combat*pred_mat_combat[,1] + coef2_combat*pred_mat_combat[,2] + coef3_combat*pred_mat_combat[,3] + coef4_combat*pred_mat_combat[,4] + coef5_combat*pred_mat_combat[,5])/(coef1_combat+coef2_combat+coef3_combat+coef4_combat+coef5_combat)
  
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = 10, nrow = nrow(testing)))
  colnames(LOSO) = c("auc_ds1","auc_ds2", "auc_ds3","auc_ds4", "auc_ds5", "auc_ds1_combat","auc_ds2_combat", "auc_ds3_combat","auc_ds4_combat", "auc_ds5_combat")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = auc(testing$status[-k], pred_prob1[-k])
    LOSO[k,2] = auc(testing$status[-k], pred_prob2[-k])
    LOSO[k,3] = auc(testing$status[-k], pred_prob3[-k])
    LOSO[k,4] = auc(testing$status[-k], pred_prob4[-k])
    LOSO[k,5] = auc(testing$status[-k], pred_prob5[-k])
    LOSO[k,6] = auc(testing$status[-k], pred_prob1_combat[-k])
    LOSO[k,7] = auc(testing$status[-k], pred_prob2_combat[-k])
    LOSO[k,8] = auc(testing$status[-k], pred_prob3_combat[-k])
    LOSO[k,9] = auc(testing$status[-k], pred_prob4_combat[-k])
    LOSO[k,10] = auc(testing$status[-k], pred_prob5_combat[-k])
  }
  LOSO[LOSO < 0.5] = 0.5
  pred_LOSO = ((LOSO[,1]-0.5)*pred_prob1 + (LOSO[,2]-0.5)*pred_prob2 + (LOSO[,3]-0.5)*pred_prob3 + (LOSO[,4]-0.5)*pred_prob4 + (LOSO[,5]-0.5)*pred_prob5)/ (LOSO[,1]+LOSO[,2]+LOSO[,3]+LOSO[,4]+LOSO[,5]-2.5)
  pred_LOSO_combat = ((LOSO[,6]-0.5)*pred_prob1_combat + (LOSO[,7]-0.5)*pred_prob2_combat + (LOSO[,8]-0.5)*pred_prob3_combat + (LOSO[,9]-0.5)*pred_prob4_combat + (LOSO[,10]-0.5)*pred_prob5_combat) / (LOSO[,6]+LOSO[,7]+LOSO[,8]+LOSO[,9]+LOSO[,10]-2.5)
  
  
  # rank methods
  rank_mat = rank_method_5trainings(testing, pred_prob1, pred_prob2, pred_prob3, pred_prob4, pred_prob5)
  rank_mat_combat = rank_method_5trainings(testing, pred_prob1_combat, pred_prob2_combat, pred_prob3_combat, pred_prob4_combat, pred_prob5_combat)
  
  
  ####  Evaluate performance
  tst_scores <- c(list(Training1 = pred_prob1, Training2 = pred_prob2, Training3 = pred_prob3, Training4 = pred_prob4, Training5 = pred_prob5, 
                       Training1_combat = pred_prob1_combat, Training2_combat = pred_prob2_combat, Training3_combat = pred_prob3_combat, Training4_combat = pred_prob4_combat, Training5_combat = pred_prob5_combat,
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
perf_df[31,] = colMeans(perf_df[1:30,])
perf_df[32,] = apply(perf_df[1:30,], 2, sd)

####  Output results
file_name <- sprintf('integration_auc_5trainings_test_%s_logRelAbun.csv', 
                     gsub('.', '', sub("/Users/lynngao/Desktop/simulation of heterogenity/TB_gene_expression/", "", sub(".rds","",command_args[6])), fixed=T))
write.csv(perf_df, paste0(command_args[7],file_name))
