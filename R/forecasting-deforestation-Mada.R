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

# Create a new environment 
virtualenv_create("r-reticulate")
# Install deforestprob
virtualenv_install("r-reticulate", "deforestprob")
# Import deforestprob
Sys.unsetenv("DISPLAY")  # Remove DISPLAY
dfp <- import("deforestprob")

# ===============================
# Download data
# ===============================

# Make data directory
dir.create("data", showWarnings=FALSE)

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
dir.create("output", showWarnings=FALSE)

# Sample points
dataset <- dfp$sample(nsamp=10000L, Seed=1234L, csize=10L,
					 var_dir="data",
					 input_forest_raster="fordefor2010.tif",
					 output_file="output/sample.txt",
					 blk_rows=1L)
# dataset <- read.table("output/sample.txt", header=TRUE, sep=",")
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
	suitability_formula=formula, data=dataset_nona,
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

# Plot traces and posteriors
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

# Install Python package in virtual environment
virtualenv_install("r-reticulate", c("patsy","statsmodels"))
patsy <- import("patsy")
sm <- import("statsmodels.api")
smf <- import("statsmodels.formula.api")

# Full model
deviance_full <- 0

# Null model
formula_null <- "I(1-fordefor2010) ~ 1"
dmat_null <- patsy$dmatrices(formula_null, data=r_to_py(dataset_nona), NA_action="drop",
											 return_type="dataframe", eval_env=-1L)
Y <- dmat_null[[1]]
X_null <- dmat_null[[2]]
mod_null <- sm$GLM(Y, X_null, family=sm$families$Binomial()).fit()
deviance_null <- mod_null$deviance

# Model with no spatial random effects
formula_nsre = paste0("I(1-fordefor2010) ~ C(sapm) + scale(altitude) + ",
											"scale(slope) + scale(dist_defor) + scale(dist_defor)*scale(dist_defor) + ",
											"scale(dist_edge) + scale(dist_road) + scale(dist_town)")
mod_nsre <- smf$glm(formula_nsre, r_to_py(dataset_nona), family=sm$families$Binomial(), eval_env=-1L)$fit()
deviance_nsre <- mod_nsre$deviance

# Model with iCAR process
deviance_icar <- mod_binomial_iCAR$Deviance

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

# ========================================================
# End
# ========================================================