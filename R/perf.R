# Libraries
require(dplyr)

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

# Accuracy_indices
accuracy_indices <- function(pred, obs) {
  
  if (identical(dim(pred),as.integer(c(2,2)))) {
    df <- pred
    n00 <- df[1,1]
    n10 <- df[2,1]
    n01 <- df[1,2]
    n11 <- df[2,2]
  } else {
    
    # Create pandas data-frame
    df <- data.frame(pred=pred, obs=obs)
    
    # Confusion matrix
    n00 <- sum((df["pred"] == 0) & (df["obs"] == 0))
    n10 <- sum((df["pred"] == 1) & (df["obs"] == 0))
    n01 <- sum((df["pred"] == 0) & (df["obs"] == 1))
    n11 <- sum((df["pred"] == 1) & (df["obs"] == 1))
  }
  
  # Accuracy indices
  N <- n11 + n10 + n00 + n01
  OA <- (n11 + n00) / N
  FOM <- n11 / (n11 + n10 + n01)
  Sensitivity <- n11 / (n11 + n01)
  Specificity <- n00 / (n00 + n10)
  TSS <- Sensitivity + Specificity - 1
  Prob_1and1 <- (n11 + n10) * (n11 + n01)
  Prob_0and0 <- (n00 + n01) * (n00 + n10)
  Expected_accuracy <- (Prob_1and1 + Prob_0and0) / (N * N)
  Kappa <- (OA - Expected_accuracy) / (1 - Expected_accuracy)
  
  r <- data.frame(OA=round(OA, 2), EA=round(Expected_accuracy, 2),
                  FOM=round(FOM, 2),
                  Sen=round(Sensitivity, 2),
                  Spe=round(Specificity, 2),
                  TSS=round(TSS, 2), K=round(Kappa, 2))
  
  return(r)
}

# Performance
performance_index <- function(data_valid, theta_pred) {
  # Model predictions for validation dataset
  data_valid$theta_pred <- theta_pred
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
    proba_thresh <- quantile(data_extract$theta_pred, 1-perc/100)  # ! must be 1-proba_defor
    data_extract$pred[data_extract$theta_pred >= proba_thresh] <- 1
    # Computing accuracy indices
    pred <- data_extract$pred
    obs <- 1-data_extract$fordefor2010
    perf <- accuracy_indices(pred,obs) %>% 
      dplyr::select(FOM,OA,EA,Spe,Sen,TSS,K)
    performance[i,2:8] <- perf
    # AUC
    pos.scores <- data_extract$theta_pred[data_extract$fordefor2010==0]
    neg.scores <- data_extract$theta_pred[data_extract$fordefor2010==1]
    performance$AUC[i] <- round(computeAUC(pos.scores,neg.scores),2)
  }
  return(performance)
}

# npix2ha
npix2ha <- function(npix, res=30) {
	return(npix*(res^2)/10000)
}

# Performance index
CV_f <- function(obs,pred) {
	RMSE <- sqrt(mean((obs-pred)^2))
	Mean <- mean(obs)
	return(RMSE/Mean)
}
R2_f <- function(obs,pred) {
	sum1 <- sum((obs-pred)^2)
	sum2 <- sum((obs-mean(obs))^2)
	return(1-sum1/sum2)
}

# Correlation as a function of scale
cor_scale <- function(df,e,square_size=33) {
	# Tidy dataset
	df <- df %>%
		dplyr::mutate(obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f),
									obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d)) %>%
		dplyr::mutate(box0=0:(nrow(df)-1)) %>%
		dplyr::filter(!(obs_f==0 & obs_d==0)) # Remove squares with no forest
	
	# Coordinates of the centers of the ~1 km boxes
	box_size0 <- 30*square_size
	ncol <- ceiling((xmax(e) - xmin(e)) / box_size0)
	nrow <- ceiling((ymax(e) - ymin(e)) / box_size0)
	y_box <- floor(df$box0/ncol) # quotient
	x_box <- df$box0 %% ncol # remainder
	y <- ymax(e)-box_size0*(y_box+0.5)
	x <- xmin(e)+box_size0*(x_box+0.5)
	# Box number from x,y coordinates and box size
	coeff <- c(1,2,5,10,15,25,50,75,100,200)
	box_size <- box_size0 * coeff # box_size in m
	nbox <- Cor <- vector()
	# Loop on box size
	for (i in 1:length(box_size)) {
		# Box size
		b <- box_size[i]
		# Number of boxes
		ncol <- ceiling((xmax(e) - xmin(e)) / b)
		nrow <- ceiling((ymax(e) - ymin(e)) / b)
		nbox[i] <- ncol * nrow
		# Box identification
		I <- floor((x - xmin(e)) / b)
		J <- floor((ymax(e) - y) / b)
		box <- J * ncol + I
		# Sum deforested areas by box
		obs_d <- df %>% mutate(box=box) %>%
			group_by(box) %>% summarise(sum_obs_d_ha=sum(obs_d_ha)) %>%
			pull(sum_obs_d_ha)
		pred_d <- df %>% mutate(box=box) %>%
			group_by(box) %>% summarise(sum_pred_d_ha=sum(pred_d_ha)) %>%
			pull(sum_pred_d_ha)
		# Correlation
		Cor[i] <- cor(obs_d,pred_d,method="pearson")
	}
	# Return
	return(data.frame(coeff=coeff,box_size=box_size,nbox=nbox,cor=Cor))
}

# End