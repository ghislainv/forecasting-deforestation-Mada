#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# Anaconda2 must be installed
# https://conda.io/docs/user-guide/install/index.html

# Libraries
require(reticulate)
require(glue)
require(curl)
require(dplyr)
require(broom)
require(ggplot2)
require(rasterVis)
require(rgdal)

# ================================
# Setup Python virtual environment
# ================================

# Force the use of miniconda2 python 2.7
#use_python("/home/ghislain/miniconda2/bin/python", required=TRUE)
#py_config()

# # Create conda virtual environment
# if (!("r-reticulate" %in% conda_list()$name)) {
# 	# Python 2.7
#   conda_create("r-reticulate","python=2.7")
# 	# Conda install with use of conda-forge
# 	conda_modules <- c("pip","gdal","numpy","scipy","statsmodels")
# 	conda_install("r-reticulate",conda_modules,forge=TRUE)
# }
# 
# # Install deforestprob from git with pip
# git_deforestprob <- "https://github.com/ghislainv/deforestprob/archive/master.zip"
# conda_install("r-reticulate",git_deforestprob,pip=TRUE,pip_ignore_installed=FALSE)
# py_discover_config("deforestprob","r-reticulate")
# 
# # Use conda env
# use_condaenv("r-reticulate", required=TRUE)

# ================================
# Import Python modules
# ================================

Sys.unsetenv("DISPLAY") # Remove DISPLAY for Python plot
dfp <- import("deforestprob") # Disregard warnings
patsy <- import("patsy")
sm <- import("statsmodels.api")
smf <- import("statsmodels.formula.api")

# ===============================
# Download data
# ===============================

# Make data directory
if (!dir.exists("data")) {
  dir.create("data")
  dir.create("data/model")
  dir.create("data/mada")
  dir.create("data/validation")
}

# Files
files <- c("fordefor2010.tif", "dist_edge.tif", "dist_defor.tif",
				 "altitude.tif", "slope.tif", "aspect.tif",
				 "dist_road.tif", "dist_town.tif", "dist_river.tif",
				 "sapm.tif", "roads.kml", "towns.kml", "rivers.kml", "sapm.kml")

# Download model data
for (i in files) {
  if (!file.exists(glue("data/model/{i}"))) {
		f <- glue("https://zenodo.org/record/259582/files/{i}")
		curl_download(f, glue("data/model/{i}"), quiet=FALSE)
	}
}

# Download validation data
if (!file.exists("data/validation/for2014.tif")) {
  f <- "http://bioscenemada.cirad.fr/FileTransfer/for2014.tif"
  curl_download(f, "data/validation/for2014.tif", quiet=FALSE)
}

# Download mada outline
if (!file.exists("data/mada/MAD_outline.shp")) {
  f <- "http://bioscenemada.cirad.fr/FileTransfer/MAD_outline.zip"
  curl_download(f, "data/mada/MAD_outline.zip", quiet=FALSE)
  unzip("data/mada/MAD_outline.zip", exdir="data/mada")
}

# ========================================================
# Sample points
# ========================================================

# Make output directory
if (!dir.exists("output")) {
  dir.create("output")
}

# Sample points
if (!file.exists("output/sample.txt")) {
  dataset <- dfp$sample(nsamp=10000L, Seed=1234L, csize=10L,
                        var_dir="data/model",
                        input_forest_raster="fordefor2010.tif",
                        output_file="output/sample.txt",
                        blk_rows=1L)
}
dataset <- read.table("output/sample.txt", header=TRUE, sep=",")
head(dataset)

# ========================================================
# hSDM model
# ========================================================

# Set number of trials to one
dataset$trials <- 1

# Remove observations with NA
dataset_nona <- dataset %>%
	dplyr::filter(complete.cases(dataset))

# Spatial cells for spatial-autocorrelation
neighborhood <- dfp$cellneigh_ctry(raster="data/model/fordefor2010.tif",
																	 vector="data/mada/MAD_outline.shp",
																	 csize=10L, rank=1L)
nneigh <- neighborhood[[1]]
adj <- neighborhood[[2]]

# Formula
formula <- paste0("I(1-fordefor2010) + trials ~ C(sapm) + scale(altitude) + scale(slope) +",
									"scale(dist_defor) + np.power(scale(dist_defor),2) + ",
									"scale(dist_edge) + ",
									"scale(dist_road) + scale(dist_town) + cell")

# Model
mod_binomial_iCAR <- dfp$model_binomial_iCAR(
	# Observations
	suitability_formula=formula, data=r_to_py(dataset_nona),
	# Spatial structure
	n_neighbors=np_array(nneigh,dtype="int32"), neighbors=np_array(adj,dtype="int32"),
	# Environment
	eval_env=-1L,
	# Chains
	burnin=1000L, mcmc=1000L, thin=1L,
	# Starting values
	beta_start=-99)

# Save model (not yet supported)
# py_save_object(mod_binomial_iCAR, "output/mod_binomial_iCAR.pkl",
# 							 pickle = "pickle")

# Summary
sink(file="output/summary_mod_binomial_iCAR.txt")
print(mod_binomial_iCAR)
sink()

# Plot traces and posteriors (try to fix warnings !!)
traces_fig <- mod_binomial_iCAR$plot(output_file="output/mcmc.pdf",
																		 plots_per_page=3L,
																		 figsize=c(9,6),
																		 dpi=100)

# ========================================================
# Resampling spatial random effects
# ========================================================

# Spatial random effects
rho <- mod_binomial_iCAR$rho

# Resample
dfp$resample_rho(rho=rho, input_raster="data/model/fordefor2010.tif",
								 output_file="output/rho.tif",
								 csize_orig=10, csize_new=1)

# Plot random effects
dfp$plot$rho("output/rho_orig.tif",output_file="output/rho_orig.png")
dfp$plot$rho("output/rho.tif",output_file="output/rho.png")

# ========================================================
# Predicting spatial probability of deforestation
# ========================================================

# Compute predictions
fig_pred <- dfp$predict_raster_binomial_iCAR(mod_binomial_iCAR, var_dir="data/model",
											  input_cell_raster="output/rho.tif",
											  input_forest_raster="data/model/fordefor2010.tif",
											  output_file="output/pred_binomial_iCAR.tif",
											  blk_rows=128L)

# Plot predictions
dfp$plot$prob("output/pred_binomial_iCAR.tif",output_file="output/pred_binomial_iCAR.png")

# ========================================================
# Predicting forest cover
# ========================================================

forest_cover <- dfp$deforest(input_raster="output/pred_binomial_iCAR.tif",
														 hectares=4000000,
														 output_file="output/forest_cover_2050.tif",
														 blk_rows=128L)

# Plot future forest cover
dfp$plot$fcc("output/forest_cover_2050.tif",output_file="output/forest_cover_2050.png")

# ========================================================
# Model comparison
# ========================================================

# Full model
deviance_full <- 0

# Null model
formula_null <- "I(1-fordefor2010) ~ 1"
dmat_null <- patsy$dmatrices(formula_null, data=r_to_py(dataset_nona), NA_action="drop",
											 return_type="dataframe", eval_env=-1L)
Y <- dmat_null[[1]]
X_null <- dmat_null[[2]]
mod_null <- sm$GLM(Y, X_null, family=sm$families$Binomial())$fit()
deviance_null <- mod_null$deviance

# Model with no spatial random effects
formula_nsre <- paste0("I(1-fordefor2010) ~ C(sapm) + scale(altitude) + ",
											"scale(slope) + scale(dist_defor) + I(scale(dist_defor)*scale(dist_defor)) + ",
											"scale(dist_edge) + scale(dist_road) + scale(dist_town)")
mod_nsre <- smf$glm(formula_nsre, r_to_py(dataset_nona), family=sm$families$Binomial(), eval_env=-1L)$fit()
deviance_nsre <- mod_nsre$deviance

# Summary nsre
sink(file="output/summary_mod_binomial_nsre.txt")
print(mod_nsre$summary())
sink()


# ================
# Accuracy indices

# 1. iCAR model
dataset_nona$pred <- 0
dataset_nona$theta_pred <- mod_binomial_iCAR$theta_pred
# Proportion of deforested pixels
nobs <- length(dataset_nona$fordefor2010)  # inferior to 20000 as there were NaN
ndefor <- sum(dataset_nona$fordefor2010==0)  # 0 for deforestation in fordefor2010
proba_defor <- ndefor/nobs  # not exactly 0.5
# Probability threshold to transform probability into binary values
proba_thresh <- quantile(dataset_nona$theta_pred, 1-proba_defor)  # ! must be (1-proba_defor)
dataset_nona$pred[dataset_nona$theta_pred >= proba_thresh] <- 1
# We check that the proportion of deforested pixel is the same for observations/predictions
ndefor_pred <- sum(dataset_nona$pred == 1)
proba_defor_pred <- ndefor_pred/nobs
proba_defor_pred == proba_defor
# Computing accuracy indices
pred <- dataset_nona$pred
obs <- 1-dataset_nona$fordefor2010
internal_validation_iCAR <- as.data.frame(dfp$accuracy_indices(pred,obs))

# 2. nsre
dataset_nona$pred <- 0
dataset_nona$theta_pred <- mod_nsre$fittedvalues
# Proportion of deforested pixels
nobs <- length(dataset_nona$fordefor2010)  # inferior to 20000 as there were NaN
ndefor <- sum(dataset_nona$fordefor2010==0)  # 0 for deforestation in fordefor2010
proba_defor <- ndefor/nobs  # not exactly 0.5
# Probability threshold to transform probability into binary values
proba_thresh <- quantile(dataset_nona$theta_pred, 1-proba_defor)  # ! must be (1-proba_defor)
dataset_nona$pred[dataset_nona$theta_pred >= proba_thresh] <- 1
# We check that the proportion of deforested pixel is the same for observations/predictions
ndefor_pred <- sum(dataset_nona$pred == 1)
proba_defor_pred <- ndefor_pred/nobs
proba_defor_pred == proba_defor
# Computing accuracy indices
pred <- dataset_nona$pred
obs <- 1-dataset_nona$fordefor2010
internal_validation_nsre <- as.data.frame(dfp$accuracy_indices(pred,obs))

# 3. With random predictions
internal_validation_rand <- matrix(data=NA,ncol=6,nrow=100)
for (i in 1:100) {
  pred <- sample(x=c(0,1),size=nobs,replace=TRUE,prob=c(1-proba_defor,proba_defor))
  internal_validation_rand[i,] <- as.numeric(dfp$accuracy_indices(pred,obs))
}
mean_rand <- apply(internal_validation_rand,2,mean)
sd_rand <- apply(internal_validation_rand,2,sd)

# 4. Combine iCAR and random
int_val <- rbind(internal_validation_iCAR,internal_validation_nsre,
                 mean_rand,sd_rand)
int_val <- round(t(int_val)*100)
colnames(int_val) <- c("iCAR","nsre","mu_rand","sd_rand")
int_val <- as.data.frame(int_val) %>%
  dplyr::mutate(index=rownames(int_val)) %>%
  dplyr::select(index,iCAR,nsre,mu_rand,sd_rand)
write.table(int_val,file="output/internal_validation.txt",
            row.names=FALSE, sep=",")
# Check here, Kappa and TSS have exact same values... might need to be corrected.

# ===============
# Model deviances

# Model with iCAR process
deviance_icar <- mod_binomial_iCAR$deviance

# Dataframe
dev <- c(deviance_null, deviance_nsre, deviance_icar, deviance_full)
mod <- data.frame(model=c("null", "nsre", "icar", "full"),
									deviance=round(dev))
perc <- 100*(1-mod$deviance/deviance_null)
mod$perc <- round(perc)
write.table(mod,file="output/deviance_model_comparison.txt",sep=",",row.names=FALSE)

# ===============
# Model deviances

# Compute predictions with logistic regression
fig_pred_nsre <- dfp$predict_raster(mod_nsre, var_dir="data/model",
																		input_forest_raster="data/model/fordefor2010.tif",
																		output_file="output/pred_nsre.tif",
																		blk_rows=128L, transform=TRUE)

# Plot predictions
dfp$plot$prob("output/pred_nsre.tif",output_file="output/pred_nsre.png")

# Future forest cover
forest_cover_nsre <- dfp$deforest(input_raster="output/pred_nsre.tif",
                                  hectares=4000000,
                                  output_file="output/forest_cover_2050_nsre.tif",
                                  blk_rows=128L)

# Plot future forest cover
dfp$plot$fcc("output/forest_cover_2050_nsre.tif",output_file="output/forest_cover_2050_nsre.png")

# Correlation between maps at ~10km
corr_10km_df <- dfp$validation_npix(r_pred="output/forest_cover_2050.tif",
                                    r_obs="output/forest_cover_2050_nsre.tif",
                                    output_file="output/npix.txt", 
                                    value_f=1,
                                    value_d=0,
                                    square_size=333L)

# Function to convert pixels to hectares
npix2ha <- function(npix,res=30) {return(npix*res*res/10000)}

# Tidy dataset
corr_10km <- corr_10km_df %>%
  dplyr::filter(!(obs_f== 0 & obs_d==0)) %>% # Remove squares with no forest
  dplyr::mutate(obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d),
                obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f))

# Check same number of forest and deforestation
f_tot_obs <- sum(corr_10km$obs_d_ha + corr_10km$obs_f_ha)
f_tot_pred <- sum(corr_10km$pred_d_ha + corr_10km$pred_f_ha)
d_tot_obs <- sum(corr_10km$obs_d_ha)
d_tot_pred <- sum(corr_10km$pred_d_ha)

# Plot
plot_pred_obs <- ggplot(corr_10km, aes(obs_d_ha,pred_d_ha)) +
  geom_point() +
  geom_smooth()
ggsave("output/plot_pred_obs.pdf", plot_pred_obs)

# Model
correlation <- round(100*cor(x=corr_10km$pred_d_ha, y=corr_10km$obs_d_ha, method=c("pearson")))
# Correlation is weak between the two maps at 10km: 62% (52% at 1km)

# Raster of differences
system(paste0("gdal_calc.py --overwrite -A output/forest_cover_2050.tif -B output/forest_cover_2050_nsre.tif ",
              "--outfile=output/diff_iCAR_nsre_2050.tif --type=Byte ",
              "--calc='255-254*(A==1)*(B==1)-255*(A==0)*(B==0)-253*(A==1)*(B==0)-253*(A==0)*(B==1)' ",
              "--co 'COMPRESS=LZW' --co 'PREDICTOR=2' ",
              "--NoDataValue=255"))

# Plot the raster of differences

# Data
r_diff <- raster("output/diff_iCAR_nsre_2050.tif")
mada.latlong <- readOGR(dsn="data/mada",layer="MAD_outline")
proj4string(mada.latlong) <- "+init=epsg:4326"
mada <- spTransform(mada.latlong,CRSobj=CRS("+init=epsg:32738"))
mada.df <- tidy(mada)

# Theme
## Setting basic theme options for plot with ggplot2
theme_base <- theme(axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    legend.position="none",
                    plot.margin=unit(c(0,0,0,0),"null"),
                    panel.spacing=unit(c(0,0,0,0),"null"),
                    plot.background=element_rect(fill="transparent"),
                    panel.background=element_rect(fill="transparent"),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.border=element_blank())

# Zooms
zoom1 <- list(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000)
zoom2 <- list(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000)
Extent.zoom1 <- paste(zoom1$xmin,zoom1$ymin,zoom1$xmax,zoom1$ymax)
Extent.zoom2 <- paste(zoom2$xmin,zoom2$ymin,zoom2$xmax,zoom2$ymax)

system(paste0("gdalwarp -overwrite -dstnodata 255 \\
          -r near -tr 30 30 \\
          -te ",Extent.zoom1," -of GTiff \\
          -co 'compress=lzw' -co 'predictor=2' \\
          output/diff_iCAR_nsre_2050.tif \\
          output/diff_zoom1.tif"))
system(paste0("gdalwarp -overwrite -dstnodata 255 \\
          -r near -tr 30 30 \\
          -te ",Extent.zoom2," -of GTiff \\
          -co 'compress=lzw' -co 'predictor=2' \\
          output/diff_iCAR_nsre_2050.tif \\
          output/diff_zoom2.tif"))

r_zoom1 <- raster("output/diff_zoom1.tif")
r_zoom2 <- raster("output/diff_zoom2.tif")

# Plot
p_diff <- gplot(r_diff,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000),
            fill="transparent",colour="black",size=0.3) +
  geom_rect(aes(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(300000,1100000),ylim=c(7165000,8685000))
ggsave("output/diff_iCAR_nsre_2050.png", p_diff, width=4.5, height=8)

# Zoom 1
zoom1 <- gplot(r_zoom1,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(346000,439000), ylim=c(7387000,7480000))
ggsave("output/diff_zoom1.png", zoom1)

# Zoom 2
zoom2 <- gplot(r_zoom2,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(793000,886000), ylim=c(7815000,7908000))
ggsave("output/diff_zoom2.png", zoom2)

# ====================================================
# Predictions vs. observations on the period 2010-2014
# ====================================================

# Compute forest-cover change between 2010 and 2014
system(paste0("gdal_translate -a_nodata 99 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' ",
              "data/model/fordefor2010.tif output/fordefor2010_.tif")) # Set nodata different from 255
system(paste0("gdal_translate -a_nodata 99 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' ",
              "data/validation/for2014.tif output/for2014_.tif"))
system(paste0("gdal_calc.py --overwrite -A output/fordefor2010_.tif -B output/for2014_.tif ",
              "--outfile=output/fcc_2010_2014_obs.tif --type=Byte ",
              "--calc='255-254*(A==1)*(B==1)-255*(A==1)*(B==255)' --co 'COMPRESS=LZW' --co 'PREDICTOR=2' ",
              "--NoDataValue=255"))

# Plot observations
dfp$plot$fcc("output/fcc_2010_2014_obs.tif",output_file="output/fcc_2010_2014_obs.png")

# Number of deforested pixels
count_d <- dfp$countpix("output/fcc_2010_2014_obs.tif", value=0L, blk_rows=128L)
count_d$area # 392601 ha

# ==============
# Comp with iCAR

# Predictions iCAR 2010-2014
fcc_iCAR_2014 <- dfp$deforest(input_raster="output/pred_binomial_iCAR.tif",
                              hectares=count_d$area,
                              output_file="output/fcc_2014_iCAR.tif",
                              blk_rows=128L)

# Correlation between maps at ~10km2
val_iCAR_df <- dfp$validation_npix(r_pred="output/fcc_2014_iCAR.tif",
                                   r_obs="output/fcc_2010_2014_obs.tif",
                                   output_file="output/npix_iCAR.txt", 
                                   value_f=1,
                                   value_d=0,
                                   square_size=333L)

# Tidy dataset
val_iCAR <- val_iCAR_df %>%
  dplyr::filter(!(obs_f== 0 & obs_d==0)) %>% # Remove squares with no forest
  dplyr::mutate(obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d),
                obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f))
write.table(val_iCAR,"output/val_iCAR.txt",row.names=FALSE,sep=",")

# Check same number of forest and deforestation
f_tot_obs <- sum(val_iCAR$obs_d_ha + val_iCAR$obs_f_ha)
f_tot_pred <- sum(val_iCAR$pred_d_ha + val_iCAR$pred_f_ha)
d_tot_obs <- sum(val_iCAR$obs_d_ha)
d_tot_pred <- sum(val_iCAR$pred_d_ha)

# Plot
pred_obs_iCAR <- ggplot(val_iCAR, aes(obs_d_ha,pred_d_ha)) +
  geom_point(alpha=0.1) +
  geom_smooth() + geom_abline(intercept=0, slope=1) +
  coord_equal(xlim=c(0,1000), ylim=c(0,1000))
ggsave("output/plot_pred_obs_iCAR.pdf", pred_obs_iCAR)

# Model
corr_iCAR <- round(100*cor(x=val_iCAR$pred_d_ha, y=val_iCAR$obs_d_ha, method=c("pearson")))

# Raster of differences
system(paste0("gdal_calc.py --overwrite -A output/fcc_2010_2014_obs.tif -B output/fcc_2014_iCAR.tif ",
              "--outfile=output/diff_obs_iCAR_2014.tif --type=Byte ",
              "--calc='255-254*(A==1)*(B==1)-255*(A==0)*(B==0)-253*(A==1)*(B==0)-253*(A==0)*(B==1)' ",
              "--co 'COMPRESS=LZW' --co 'PREDICTOR=2' ",
              "--NoDataValue=255"))

system(paste0("gdalwarp -overwrite -dstnodata 255 \\
          -r near -tr 30 30 \\
          -te ",Extent.zoom1," -of GTiff \\
          -co 'compress=lzw' -co 'predictor=2' \\
          output/diff_obs_iCAR_2014.tif \\
          output/diff_zoom1_iCAR.tif"))
system(paste0("gdalwarp -overwrite -dstnodata 255 \\
          -r near -tr 30 30 \\
          -te ",Extent.zoom2," -of GTiff \\
          -co 'compress=lzw' -co 'predictor=2' \\
          output/diff_obs_iCAR_2014.tif \\
          output/diff_zoom2_iCAR.tif"))

r_zoom1_iCAR <- raster("output/diff_zoom1_iCAR.tif")
r_zoom2_iCAR <- raster("output/diff_zoom2_iCAR.tif")

# Plot the raster of differences

# Data
r_diff_iCAR_2014 <- raster("output/diff_obs_iCAR_2014.tif")

# Plot
p_diff_iCAR_2014 <- gplot(r_diff_iCAR_2014,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000),
            fill="transparent",colour="black",size=0.3) +
  geom_rect(aes(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(300000,1100000),ylim=c(7165000,8685000))
ggsave("output/diff_obs_iCAR_2014.png", p_diff_iCAR_2014, width=4.5, height=8)

# Zoom 1
zoom1_iCAR <- gplot(r_zoom1_iCAR,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=346000,xmax=439000,ymin=7387000,ymax=7480000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(346000,439000), ylim=c(7387000,7480000))
ggsave("output/diff_zoom1_iCAR.png", zoom1_iCAR)

# Zoom 2
zoom2_iCAR <- gplot(r_zoom2_iCAR,maxpixels=10e5) +
  geom_polygon(data=mada.df, aes(x=long, y=lat, group=id), colour="black", fill="white", size=0.3) +
  geom_raster(aes(fill=factor(value))) +
  scale_fill_manual(values=c("red","forestgreen","blue"), na.value="transparent") +
  geom_rect(aes(xmin=793000,xmax=886000,ymin=7815000,ymax=7908000),
            fill="transparent",colour="black",size=0.3) +
  theme_bw() + theme_base +
  coord_equal(xlim=c(793000,886000), ylim=c(7815000,7908000))
ggsave("output/diff_zoom2_iCAR.png", zoom2_iCAR)

# ==============
# Comp with nsre

# Predictions nsre 2010-2014
fcc_nsre_2014 <- dfp$deforest(input_raster="output/pred_nsre.tif",
                              hectares=count_d$area,
                              output_file="output/fcc_2014_nsre.tif",
                              blk_rows=128L)

# Correlation between maps at ~10km2
val_nsre_df <- dfp$validation_npix(r_pred="output/fcc_2014_nsre.tif",
                                   r_obs="output/fcc_2010_2014_obs.tif",
                                   output_file="output/npix_nsre.txt", 
                                   value_f=1,
                                   value_d=0,
                                   square_size=333L)

# Tidy dataset
val_nsre <- val_nsre_df %>%
  dplyr::filter(!(obs_f== 0 & obs_d==0)) %>% # Remove squares with no forest
  dplyr::mutate(obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d),
                obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f))
write.table(val_nsre,"output/val_nsre.txt",row.names=FALSE,sep=",")

# Plot
pred_obs_nsre <- ggplot(val_nsre, aes(obs_d_ha,pred_d_ha)) +
  geom_point(alpha=0.1) +
  geom_smooth() + geom_abline(intercept=0, slope=1) + 
  coord_equal(xlim=c(0,1000), ylim=c(0,1000))
ggsave("output/plot_pred_obs_nsre.pdf", pred_obs_nsre)

# Model
corr_nsre <- round(100*cor(x=val_nsre$pred_d_ha, y=val_nsre$obs_d_ha, method=c("pearson")))

# Backup correlation estimates
corr_df <- data.frame(model=c("nsre","iCAR"),
                     corr=c(corr_nsre,corr_iCAR))
write.table(corr_df,"output/corr_df.txt",row.names=FALSE,sep=",")

# ========================================================
# End
# ========================================================