#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# Libraries
require(reticulate)
require(glue)
require(curl)
require(dplyr)
require(ggplot2)

# ================================
# Setup Python virtual environment
# ================================

# Set directory for virtual env
Sys.setenv(WORKON_HOME=getwd())

# Create virtual env if necessary
if (!dir.exists("r-reticulate")) {
  # Create
  virtualenv_create("r-reticulate")
  # Install Python modules
  virtualenv_install("r-reticulate",
                     c("patsy","statsmodels"))
  # virtualenv_remove("r-reticulate","deforestprob")
  virtualenv_install("r-reticulate","deforestprob",
                     ignore_installed=TRUE)
}

# Use virtual env
use_virtualenv("r-reticulate")

# Import Python modules
Sys.unsetenv("DISPLAY") # Remove DISPLAY for Python plot
dfp <- import("deforestprob")
patsy <- import("patsy")
sm <- import("statsmodels.api")
smf <- import("statsmodels.formula.api")

# ===============================
# Download data
# ===============================

# Make data directory
if (!dir.exists("data")) {
  dir.create("data")
}

# Files
files <- c("fordefor2010.tif", "dist_edge.tif", "dist_defor.tif",
				 "altitude.tif", "slope.tif", "aspect.tif",
				 "dist_road.tif", "dist_town.tif", "dist_river.tif",
				 "sapm.tif", "roads.kml", "towns.kml", "rivers.kml", "sapm.kml")

# Download
for (i in files) {
  if (!file.exists(glue("data/{i}"))) {
		f <- glue("https://zenodo.org/record/259582/files/{i}")
		curl_download(f, glue("data/{i}"), quiet=FALSE)
	}
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
                        var_dir="data",
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
neighborhood <- dfp$cellneigh(raster="data/fordefor2010.tif", csize=10L, rank=1L)
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
dfp$resample_rho(rho=rho, input_raster="data/fordefor2010.tif",
								 output_file="output/rho.tif",
								 csize_orig=10, csize_new=1)

# Plot random effects
dfp$plot$rho("output/rho_orig.tif",output_file="output/rho_orig.png")
dfp$plot$rho("output/rho.tif",output_file="output/rho.png")

# ========================================================
# Predicting spatial probability of deforestation
# ========================================================

# Compute predictions
fig_pred <- dfp$predict_raster_binomial_iCAR(mod_binomial_iCAR, var_dir="data",
											  input_cell_raster="output/rho.tif",
											  input_forest_raster="data/fordefor2010.tif",
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
formula_nsre = paste0("I(1-fordefor2010) ~ C(sapm) + scale(altitude) + ",
											"scale(slope) + scale(dist_defor) + scale(dist_defor)*scale(dist_defor) + ",
											"scale(dist_edge) + scale(dist_road) + scale(dist_town)")
mod_nsre <- smf$glm(formula_nsre, r_to_py(dataset_nona), family=sm$families$Binomial(), eval_env=-1L)$fit()
deviance_nsre <- mod_nsre$deviance

# Model with iCAR process
deviance_icar <- mod_binomial_iCAR$deviance

# Dataframe
dev <- c(deviance_null, deviance_nsre, deviance_icar, deviance_full)
mod <- data.frame(model=c("null", "nsre", "icar", "full"),
									deviance=round(dev))
perc <- 100*(1-mod$deviance/deviance_null)
mod$perc <- round(perc)
write.table(mod,file="output/deviance_model_comparison.txt",sep=",",row.names=FALSE)

# Compute predictions with logistic regression
fig_pred_nsre <- dfp$predict_raster(mod_nsre, var_dir="data",
																		input_forest_raster="data/fordefor2010.tif",
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

# ====================================================
# Predictions vs. observations on the period 2010-2014
# ====================================================

# Download
if (!file.exists("output/for2014.tif")) {
  f <- "http://bioscenemada.cirad.fr/FileTransfer/for2014.tif"
  curl_download(f, "output/for2014.tif", quiet=FALSE)
}

# Compute forest-cover change between 2010 and 2014
system(paste0("gdal_translate -a_nodata 99 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' ",
                          "data/fordefor2010.tif output/fordefor2010_.tif")) # Set nodata different from 255
system(paste0("gdal_translate -a_nodata 99 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' ",
                          "output/for2014.tif output/for2014_.tif"))
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

# Correlation between maps at ~20km2
val_iCAR_df <- dfp$validation_npix(r_pred="output/fcc_2014_iCAR.tif",
                                   r_obs="output/fcc_2010_2014_obs.tif",
                                   output_file="output/npix_iCAR.txt", 
                                   value_f=1,
                                   value_d=0,
                                   square_size=666L)

# Tidy dataset
val_iCAR <- val_iCAR_df %>%
  dplyr::filter(!(obs_f== 0 & obs_d==0)) %>% # Remove squares with no forest
  dplyr::mutate(obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d),
                obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f))

# Plot
pred_obs_iCAR <- ggplot(val_iCAR, aes(obs_d_ha,pred_d_ha)) +
  geom_point() +
  geom_smooth()
ggsave("output/plot_pred_obs_iCAR.pdf", pred_obs_iCAR)

# Model
corr_iCAR <- round(100*cor(x=val_iCAR$pred_d_ha, y=val_iCAR$obs_d_ha, method=c("pearson")))

# ==============
# Comp with nsre

# Predictions nsre 2010-2014
fcc_nsre_2014 <- dfp$deforest(input_raster="output/pred_nsre.tif",
                              hectares=count_d$area,
                              output_file="output/fcc_2014_nsre.tif",
                              blk_rows=128L)

# Correlation between maps at ~10km
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

# Plot
pred_obs_nsre <- ggplot(val_nsre, aes(obs_d_ha,pred_d_ha)) +
  geom_point() +
  geom_smooth()
ggsave("output/plot_pred_obs_nsre.pdf", pred_obs_nsre)

# Model
corr_nsre <- round(100*cor(x=val_nsre$pred_d_ha, y=val_nsre$obs_d_ha, method=c("pearson")))

# Backup correlation estimates
corr_df <- data.frame(model=c("nsre","iCAR"),
                     corr=c(corr_nsre,corr_iCAR))
write.table(corr_df,"output/corr_df",row.names=FALSE,sep=",")

# ========================================================
# End
# ========================================================