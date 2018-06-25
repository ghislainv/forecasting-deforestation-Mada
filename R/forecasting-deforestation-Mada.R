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
dfp <- import("deforestprob")
np <- import("numpy")

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
if (!file.exists(glue("data/{i}"))) {
	for (i in files) {
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

patsy <- import("patsy")
patsy$dmatrices("fordefor2010~sapm", dataset_nona, "drop")

# Model
mod_binomial_iCAR <- dfp$model_binomial_iCAR(
	# Observations
	suitability_formula=formula, data=dataset_nona,
	# Spatial structure
	n_neighbors=np_array(nneigh,dtype="int32"), neighbors=np_array(adj,dtype="int32"),
	# Chains
	burnin=1000L, mcmc=1000L, thin=1L,
	# Starting values
	beta_start=r_to_py(-99))

