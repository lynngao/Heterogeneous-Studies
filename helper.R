### ComBat normalization
combat_norm <- function(ds1, ds2){
  merged = rbind(ds1, ds2)
  data = t(merged[,-ncol(merged)])
  batch = c(rep(0,nrow(ds1)), rep(1,nrow(ds2)))
  data_combat = as.data.frame(t(ComBat(data, batch, ref.batch = 0)))
  ds2_combat = data_combat[(nrow(ds1)+1):nrow(data_combat),]
  ds2_combat$status = ds2$status
  return(ds2_combat)  
}

### Train ML model and predict
ml_model <- function(training, method) {
  if (method == "rf"){
    model <-  train(
      status ~ ., data = training, method = "ranger", num.trees = 1000, metric = "ROC",
      trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
    )
  }
  if (method == "logit"){
    parameters <- seq(0,1,0.001)
    model<- train(
      status ~ ., data = training, method = "glmnet", family="binomial",metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
      trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
    )
  }
  return(model)
}

### five rank methods: mean, geometric mean, Stuart, RRA, Bayesian rank analysis
library(RobustRankAggreg)
library(BayesRankAnalysis)
rank_method <- function(test_expr_norm, pred_prob1, pred_prob2){
  rank_mat = as.data.frame(matrix(NA, nrow = nrow(test_expr_norm), ncol = 11))
  rownames(rank_mat) = rownames(test_expr_norm)
  colnames(rank_mat) = c("prob1","prob2","rank1","rank2","rank_ratio1","rank_ratio2","mean","geo_mean","Q","rho","Bayesian")
  rank_mat$prob1 = pred_prob1
  rank_mat$prob2 = pred_prob2
  rank_mat$rank1 = rank(-rank_mat$prob1, ties.method="min")
  rank_mat$rank2 = rank(-rank_mat$prob2, ties.method="min")
  rank_mat$rank_ratio1 = rank_mat$rank1/nrow(rank_mat)
  rank_mat$rank_ratio2 = rank_mat$rank2/nrow(rank_mat)
  ##mean and geometric mean
  geo_mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2")], N = nrow(rank_mat), method = "geom.mean")
  mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2")], N = nrow(rank_mat), method = "mean")
  ##Stuart and RRA
  Q = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2")], N = nrow(rank_mat), method = "stuart")
  rho = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2")], N = nrow(rank_mat))
  for (l in 1:nrow(rank_mat)){
    rank_mat$geo_mean[l] = geo_mean$Score[geo_mean$Name == rownames(rank_mat)[l]]
    rank_mat$mean[l] = mean$Score[mean$Name == rownames(rank_mat)[l]]
    rank_mat$Q[l] = Q$Score[Q$Name == rownames(rank_mat)[l]]
    rank_mat$rho[l] = rho$Score[rho$Name == rownames(rank_mat)[l]]
  }
  ##Bayesian rank
  M = 2  ## number of rankers
  N = nrow(rank_mat)  ## number of ranked items
  fullrank.real = as.matrix(cbind(rank_mat$rank1,rank_mat$rank2))  ## observed ranking lists
  pair.comp = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
  for(m in 1:M){
    pair.comp[,,m] = FullRankToPairComp(fullrank.real[,m] )
  }
  iter.max = 1000   ## Gibbs sampler total iterations
  iter.burn = 200   ## Gibbs sampler burn-in iterations
  print.opt = 100  ## print a message every print.opt steps
  BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp, X.mat = matrix(NA, nrow =dim(pair.comp)[1], ncol = 0), 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
  BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
  rank_mat$Bayesian = BAR.fit$agg.rank
  ##BIRRA
  #rank_mat$BIRRA = BIRRA(fullrank.real)
  return(rank_mat)
}

rank_method_3trainings <- function(test_expr_norm, pred_prob1, pred_prob2, pred_prob3){
  rank_mat = as.data.frame(matrix(NA, nrow = nrow(test_expr_norm), ncol = 14))
  rownames(rank_mat) = rownames(test_expr_norm)
  colnames(rank_mat) = c("prob1","prob2","prob3",
                         "rank1","rank2","rank3",
                         "rank_ratio1","rank_ratio2","rank_ratio3",
                         "mean","geo_mean","Q","rho","Bayesian")
  rank_mat$prob1 = pred_prob1
  rank_mat$prob2 = pred_prob2
  rank_mat$prob3 = pred_prob3
  rank_mat$rank1 = rank(-rank_mat$prob1, ties.method="min")
  rank_mat$rank2 = rank(-rank_mat$prob2, ties.method="min")
  rank_mat$rank3 = rank(-rank_mat$prob3, ties.method="min")
  rank_mat$rank_ratio1 = rank_mat$rank1/nrow(rank_mat)
  rank_mat$rank_ratio2 = rank_mat$rank2/nrow(rank_mat)
  rank_mat$rank_ratio3 = rank_mat$rank3/nrow(rank_mat)
  ##mean and geometric mean
  geo_mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3")], N = nrow(rank_mat), method = "geom.mean")
  mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3")], N = nrow(rank_mat), method = "mean")
  ##Stuart and RRA
  Q = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3")], N = nrow(rank_mat), method = "stuart")
  rho = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3")], N = nrow(rank_mat))
  for (l in 1:nrow(rank_mat)){
    rank_mat$geo_mean[l] = geo_mean$Score[geo_mean$Name == rownames(rank_mat)[l]]
    rank_mat$mean[l] = mean$Score[mean$Name == rownames(rank_mat)[l]]
    rank_mat$Q[l] = Q$Score[Q$Name == rownames(rank_mat)[l]]
    rank_mat$rho[l] = rho$Score[rho$Name == rownames(rank_mat)[l]]
  }
  ##Bayesian rank
  M = 3  ## number of rankers
  N = nrow(rank_mat)  ## number of ranked items
  fullrank.real = as.matrix(cbind(rank_mat$rank1,rank_mat$rank2,rank_mat$rank3))  ## observed ranking lists
  pair.comp = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
  for(m in 1:M){
    pair.comp[,,m] = FullRankToPairComp(fullrank.real[,m] )
  }
  iter.max = 1000   ## Gibbs sampler total iterations
  iter.burn = 200   ## Gibbs sampler burn-in iterations
  print.opt = 100  ## print a message every print.opt steps
  BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp, X.mat = matrix(NA, nrow =dim(pair.comp)[1], ncol = 0), 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
  BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
  rank_mat$Bayesian = BAR.fit$agg.rank
  ##BIRRA
  #rank_mat$BIRRA = BIRRA(fullrank.real)
  return(rank_mat)
}

rank_method_4trainings <- function(test_expr_norm, pred_prob1, pred_prob2, pred_prob3, pred_prob4){
  rank_mat = as.data.frame(matrix(NA, nrow = nrow(test_expr_norm), ncol = 17))
  rownames(rank_mat) = rownames(test_expr_norm)
  colnames(rank_mat) = c("prob1","prob2","prob3","prob4",
                         "rank1","rank2","rank3","rank4",
                         "rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4",
                         "mean","geo_mean","Q","rho","Bayesian")
  rank_mat$prob1 = pred_prob1
  rank_mat$prob2 = pred_prob2
  rank_mat$prob3 = pred_prob3
  rank_mat$prob4 = pred_prob4
  rank_mat$rank1 = rank(-rank_mat$prob1, ties.method="min")
  rank_mat$rank2 = rank(-rank_mat$prob2, ties.method="min")
  rank_mat$rank3 = rank(-rank_mat$prob3, ties.method="min")
  rank_mat$rank4 = rank(-rank_mat$prob4, ties.method="min")
  rank_mat$rank_ratio1 = rank_mat$rank1/nrow(rank_mat)
  rank_mat$rank_ratio2 = rank_mat$rank2/nrow(rank_mat)
  rank_mat$rank_ratio3 = rank_mat$rank3/nrow(rank_mat)
  rank_mat$rank_ratio4 = rank_mat$rank4/nrow(rank_mat)
  ##mean and geometric mean
  geo_mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4")], N = nrow(rank_mat), method = "geom.mean")
  mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4")], N = nrow(rank_mat), method = "mean")
  ##Stuart and RRA
  Q = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4")], N = nrow(rank_mat), method = "stuart")
  rho = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4")], N = nrow(rank_mat))
  for (l in 1:nrow(rank_mat)){
    rank_mat$geo_mean[l] = geo_mean$Score[geo_mean$Name == rownames(rank_mat)[l]]
    rank_mat$mean[l] = mean$Score[mean$Name == rownames(rank_mat)[l]]
    rank_mat$Q[l] = Q$Score[Q$Name == rownames(rank_mat)[l]]
    rank_mat$rho[l] = rho$Score[rho$Name == rownames(rank_mat)[l]]
  }
  ##Bayesian rank
  M = 4  ## number of rankers
  N = nrow(rank_mat)  ## number of ranked items
  fullrank.real = as.matrix(cbind(rank_mat$rank1,rank_mat$rank2,rank_mat$rank3,rank_mat$rank4))  ## observed ranking lists
  pair.comp = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
  for(m in 1:M){
    pair.comp[,,m] = FullRankToPairComp(fullrank.real[,m] )
  }
  iter.max = 1000   ## Gibbs sampler total iterations
  iter.burn = 200   ## Gibbs sampler burn-in iterations
  print.opt = 100  ## print a message every print.opt steps
  BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp, X.mat = matrix(NA, nrow =dim(pair.comp)[1], ncol = 0), 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
  BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
  rank_mat$Bayesian = BAR.fit$agg.rank
  ##BIRRA
  #rank_mat$BIRRA = BIRRA(fullrank.real)
  return(rank_mat)
}

rank_method_5trainings <- function(test_expr_norm, pred_prob1, pred_prob2, pred_prob3, pred_prob4, pred_prob5){
  rank_mat = as.data.frame(matrix(NA, nrow = nrow(test_expr_norm), ncol = 20))
  rownames(rank_mat) = rownames(test_expr_norm)
  colnames(rank_mat) = c("prob1","prob2","prob3","prob4","prob5",
                         "rank1","rank2","rank3","rank4","rank5",
                         "rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4","rank_ratio5",
                         "mean","geo_mean","Q","rho","Bayesian")
  rank_mat$prob1 = pred_prob1
  rank_mat$prob2 = pred_prob2
  rank_mat$prob3 = pred_prob3
  rank_mat$prob4 = pred_prob4
  rank_mat$prob5 = pred_prob5
  rank_mat$rank1 = rank(-rank_mat$prob1, ties.method="min")
  rank_mat$rank2 = rank(-rank_mat$prob2, ties.method="min")
  rank_mat$rank3 = rank(-rank_mat$prob3, ties.method="min")
  rank_mat$rank4 = rank(-rank_mat$prob4, ties.method="min")
  rank_mat$rank5 = rank(-rank_mat$prob5, ties.method="min")
  rank_mat$rank_ratio1 = rank_mat$rank1/nrow(rank_mat)
  rank_mat$rank_ratio2 = rank_mat$rank2/nrow(rank_mat)
  rank_mat$rank_ratio3 = rank_mat$rank3/nrow(rank_mat)
  rank_mat$rank_ratio4 = rank_mat$rank4/nrow(rank_mat)
  rank_mat$rank_ratio5 = rank_mat$rank5/nrow(rank_mat)
  ##mean and geometric mean
  geo_mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4","rank_ratio5")], N = nrow(rank_mat), method = "geom.mean")
  mean = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4","rank_ratio5")], N = nrow(rank_mat), method = "mean")
  ##Stuart and RRA
  Q = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4","rank_ratio5")], N = nrow(rank_mat), method = "stuart")
  rho = aggregateRanks(rmat = rank_mat[,c("rank_ratio1","rank_ratio2","rank_ratio3","rank_ratio4","rank_ratio5")], N = nrow(rank_mat))
  for (l in 1:nrow(rank_mat)){
    rank_mat$geo_mean[l] = geo_mean$Score[geo_mean$Name == rownames(rank_mat)[l]]
    rank_mat$mean[l] = mean$Score[mean$Name == rownames(rank_mat)[l]]
    rank_mat$Q[l] = Q$Score[Q$Name == rownames(rank_mat)[l]]
    rank_mat$rho[l] = rho$Score[rho$Name == rownames(rank_mat)[l]]
  }
  ##Bayesian rank
  M = 5  ## number of rankers
  N = nrow(rank_mat)  ## number of ranked items
  fullrank.real = as.matrix(cbind(rank_mat$rank1,rank_mat$rank2,rank_mat$rank3,rank_mat$rank4,rank_mat$rank5))  ## observed ranking lists
  pair.comp = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
  for(m in 1:M){
    pair.comp[,,m] = FullRankToPairComp(fullrank.real[,m] )
  }
  iter.max = 1000   ## Gibbs sampler total iterations
  iter.burn = 200   ## Gibbs sampler burn-in iterations
  print.opt = 100  ## print a message every print.opt steps
  BAR.fit = BayesRankCovSimp(pair.comp.ten = pair.comp, X.mat = matrix(NA, nrow =dim(pair.comp)[1], ncol = 0), 
                             tau2.alpha = 1^2, nu.alpha = 3,
                             tau2.beta = 10^2, nu.beta = 3,
                             iter.max = iter.max, print.opt = print.opt)
  BAR.fit$agg.rank = apply(BAR.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
  rank_mat$Bayesian = BAR.fit$agg.rank
  ##BIRRA
  #rank_mat$BIRRA = BIRRA(fullrank.real)
  return(rank_mat)
}


####  Helpers related to inverse gamma distribution
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}

ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}

####  Simulate batch effects
simBatch <- function(curr_train_expr, N_sample_size, batches_ind, batch, hyper_pars){
  n_batches <- N_sample_size# number of samples in each batch
  n_genes <- ncol(curr_train_expr)-1
  
  ## Organize hyper batch parameters
  batch_par <- list()
  for(i in 1:length(n_batches)){
    batch_par[[i]] <- sapply(hyper_pars, function(item){item[i]}) # mean, sd of gaussian; alpha, beta for InvGamma
  }
  
  ## Simulate batch parameters from hyper-pars
  gamma <- delta2 <- list()
  for(i in 1:length(n_batches)){
    gamma[[i]] <- rnorm(n_genes, mean=batch_par[[i]]["hyper_mu"], sd=batch_par[[i]]["hyper_sd"])
    delta2[[i]] <- rinvgamma(n_genes, shape=batch_par[[i]]["hyper_alpha"], scale=batch_par[[i]]["hyper_beta"])
  }
  
  ## Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~curr_train_expr$status, data=data.frame(Condition=curr_train_expr$status))
  beta <- solve(t(X) %*% X) %*% t(X) %*% as.matrix(curr_train_expr[,-ncol(curr_train_expr)])
  resid <- curr_train_expr[,-ncol(curr_train_expr)] - X %*% beta
  #range(apply(resid,1,mean)); range(apply(resid,1,var))
  
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- matrix(NA, nrow=ncol(curr_train_expr)-1, ncol=nrow(curr_train_expr), dimnames=dimnames(t(curr_train_expr[,-ncol(curr_train_expr)])))
  for(j in 1:length(n_batches)){
    curr_resid <- resid[as.character(batches_ind[[j]]),]
    spikein_var <- lapply(1:n_batches[j], function(row_ind){curr_resid[row_ind, ] * sqrt(delta2[[j]])})
    resid_varbatch[,as.character(batches_ind[[j]])] <- t(as.matrix(do.call(rbind,spikein_var)))
  }
  resid_varbatch = t(resid_varbatch)
  #sapply(1:5, function(k){mean(apply(resid[, batches_ind[[k]]],1,var))})
  #sapply(1:5, function(k){mean(apply(resid_varbatch[, batches_ind[[k]]],1,var))})
  
  # construct mean batch parameter design matrix using gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=factor(batch)))
  gamma_vec <- do.call(rbind, gamma)  #apply(gamma_vec,1,mean)
  
  # new data with added batch effect
  new_dat <- as.data.frame(cbind(X, X_batch) %*% rbind(beta, gamma_vec) + resid_varbatch)
  new_dat$status = curr_train_expr$status
  rownames(new_dat) <- rownames(curr_train_expr)
  colnames(new_dat) <- colnames(curr_train_expr)
  
  res <- list(new_dat=new_dat, batch_par=batch_par)
  return(res)
}


CS_zmatrix <- function(n_batch, training, perf_name, method){
  zmat <- matrix(0, nrow=n_batch, ncol=n_batch)
  for(i in 1:n_batch){
    for(j in 1:n_batch){
      if(i!=j){
        if (method == "rf"){
          model <- train(
            status ~ ., data = training[[i]], method = "ranger", num.trees = 1000, metric = "ROC",
            trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
          )
        }
        if (method == "logit"){
          parameters <- seq(0,0.5,0.001)
          model<- train(
            status ~ ., data = training[[i]], method = "glmnet", family="binomial",metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
            trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
          )
        }
        tmp_pred <- predict(model, training[[j]], type="prob")$case
        if(perf_name=="mxe"){tmp_pred <- pmax(pmin(tmp_pred, 1 - 1e-15), 1e-15)}
        # avoid Inf in computing cross-entropy loss
        test_label = ifelse(training[[j]]$status == "ctr", 0, 1)
        rocr_pred <- prediction(tmp_pred, test_label)
        perf_tmp <- performance(rocr_pred, perf_name)  # mean cross entropy
        zmat[i,j] <- as.numeric(perf_tmp@y.values)
      }
    }
  }
  return(zmat)
}
# calculate weights for CS-Avg
CS_weight <- function(cs_zmat){
  z_seq <- sqrt(rowSums(cs_zmat) / (nrow(cs_zmat)-1))
  weights_seq <- abs(z_seq - max(z_seq))
  if(sum(weights_seq)==0){weights_seq <- rep(1, length(z_seq))}
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}

Reg_SSL_pred <- function(n_batch, training, method){
  SSL_pred_lst <- SSL_coef_lst <- list()
  for(k in 1:n_batch){
    tmp <- lapply(1:n_batch, function(batch_id){
      if (method == "rf"){
        model <- train(
          status ~ ., data = training[[batch_id]], method = "ranger", num.trees = 1000, metric = "ROC",
          trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
        )
      }
      if (method == "logit"){
        parameters <- seq(0,0.5,0.001)
        model <- train(
          status ~ ., data = training[[batch_id]], method = "glmnet", family="binomial",metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
          trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
        )
      }
      pred <- predict(model, training[[k]], type="prob")
      return(pred$case)
    })
    names(tmp) <- paste0("Simulation", 1:n_batch)
    SSL_pred_lst[[k]] <- do.call(cbind, tmp)
    test_label = ifelse(training[[k]]$status == "ctr", 0, 1)
    coef_k <- nnls(A=SSL_pred_lst[[k]], b=test_label)$x  
    SSL_coef_lst[[k]] <- coef_k
  }
  return(list(pred=SSL_pred_lst, coef=SSL_coef_lst))
}

Reg_a_weight <- function(coef_mat, n_seq){
  weights_seq <- rep(0, length(n_seq))
  for(i in 1:length(n_seq)){
    weights_seq <- weights_seq + coef_mat[i, ] * n_seq[i]
  }
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}
