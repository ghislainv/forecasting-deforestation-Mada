# AUC
computeAUC <- function(pos.scores, neg.scores, n_sample=100000) {
  # Args:
  #   pos.scores: scores of positive observations
  #   neg.scores: scores of negative observations
  #   n_samples : number of samples to approximate AUC
  
  pos.sample <- sample(pos.scores, n_sample, replace=TRUE)
  neg.sample <- sample(neg.scores, n_sample, replace=TRUE)
  AUC <- mean(1.0*(pos.sample > neg.sample) + 0.5*(pos.sample==neg.sample))
  return(AUC)
}

# Performance
performance_index <- function(data_valid, model="nochange") {
  # Model predictions for validation dataset
  if (model=="nochange") {
    data_valid$theta_pred <- 0
  } else if (model=="null") {
    data_valid$theta_pred <- runif(nrow(data_valid))
  } else {
    data_valid$theta_pred <- model$predict(data_valid)
  }
  # Number of observations
  nforest <- sum(data_valid$fordefor2010==1)  # 1 for forest in fordefor2010
  ndefor <- sum(data_valid$fordefor2010==0)
  which_forest <- which(data_valid$fordefor2010==1)							
  which_defor <- which(data_valid$fordefor2010==0)
  # Performance at 1%, 10%, 25%, 50% change
  performance <- data.frame(perc=c(1,5,10,25,50),FOM=NA,OA=NA,EA=NA,
                            Spe=NA,Sen=NA,TSS=NA,K=NA,AUC=NA)
  # Loop on prevalence
  for (i in 1:length(performance$perc)) {
    perc <- performance$perc[i]
    ndefor_samp <- min(round(nforest*(perc/100)/(1-perc/100)),ndefor)
    samp_defor <- sample(which_defor,size=ndefor_samp,replace=FALSE)
    data_extract <- data_valid[c(which_forest,samp_defor),]
    data_extract$pred <- 0
    # Probability threshold to transform probability into binary values
    if (model=="nochange") {
      proba_thresh <- 1
    } else {
      proba_thresh <- quantile(data_extract$theta_pred, 1-perc/100)  # ! must be 1-proba_defor
    }
    data_extract$pred[data_extract$theta_pred >= proba_thresh] <- 1
    # Computing accuracy indices
    pred <- data_extract$pred
    obs <- 1-data_extract$fordefor2010
    perf <- as.data.frame(far$accuracy_indices(pred,obs)) %>% 
      dplyr::select(FOM,OA,EA,Spe,Sen,TSS,K)
    performance[i,2:8] <- perf
    # AUC
    pos.scores <- data_extract$theta_pred[data_extract$fordefor2010==0]
    neg.scores <- data_extract$theta_pred[data_extract$fordefor2010==1]
    performance$AUC[i] <- round(computeAUC(pos.scores,neg.scores),2)
  }
  return(performance)
}

# End