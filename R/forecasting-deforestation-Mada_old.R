#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# Environmental variables
Sys.unsetenv("DISPLAY") # Remove DISPLAY for Python plot

# Libraries
require(reticulate)
require(glue)
require(curl)
require(dplyr)
require(broom)
require(ggplot2)
require(rasterVis)
require(rgdal)

# Source R plotting functions
source("R/dfp_plot.R")

# ================================
# Setup Python virtual environment
# ================================

# Note:
# Linux: install python, gdal and virtualenv
# Windows: install miniconda: https://conda.io/docs/user-guide/install/index.html

# forestatrisk develoment version
git_forestatrisk <- "https://github.com/ghislainv/forestatrisk/archive/master.zip"

# Linux/Windows systems
if (R.Version()=="linux-gnu") {
  if (!dir.exists("far_venv")) {
    # Local virtualenv
    Sys.setenv(WORKON_HOME=getwd())
    # Create virtualenv
    virtualenv_create("far_venv")
    # Install forestatrisk package
    virtualenv_install("far_venv",git_forestatrisk,ignore_installed=TRUE)
    virtualenv_install("far_venv","statsmodels")
  }
  # Activate this virtualenv
  use_virtualenv("far_venv", required=TRUE)
} else {
  if (!("r-reticulate" %in% conda_list()$name)) {
    # Create conda virtual environment
    conda_create("r-reticulate","python=3.7")
    # Conda install with use of conda-forge
    conda_modules <- c("pip","gdal","numpy","scipy","statsmodels")
    conda_install("r-reticulate",conda_modules,forge=TRUE)
    conda_install("r-reticulate",git_forestatrisk,pip=TRUE,pip_ignore_installed=FALSE)
  }
  # Activate this virtualenv
  use_condaenv("r-reticulate", required=TRUE)
}

# Check python version and virtualenv
py_config()

# ================================
# Import Python modules
# ================================

far <- import("forestatrisk")
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
if (!file.exists("data/validation/for2017.tif")) {
  f <- "http://bioscenemada.cirad.fr/FileTransfer/for2014.tif"
  curl_download(f, "data/validation/for2014.tif", quiet=FALSE)
}

# Download mada outline
if (!file.exists("data/mada/mada38s.shp")) {
  f <- "http://bioscenemada.cirad.fr/FileTransfer/mada38s.zip"
  curl_download(f, "data/mada/mada38s.zip", quiet=FALSE)
  unzip("data/mada/mada38s.zip", exdir="data/mada")
}

# ========================================================
# Sample points
# ========================================================

# Make output directory
if (!dir.exists("output")) {
  dir.create("output")
}

# Training data-set
if (!file.exists("output/sample.txt")) {
  samp <- far$sample(nsamp=10000L, Seed=1234L, csize=10L,
                     var_dir="data/model",
                     input_forest_raster="fordefor2010.tif",
                     output_file="output/sample.txt",
                     blk_rows=1L)
}
samp <- read.table("output/sample.txt", header=TRUE, sep=",")
set.seed(1234)
train <- sample(1:20000, size=10000, replace=FALSE)
data_train <- samp[train,] %>% dplyr::filter(complete.cases(.))
data_valid <- samp[-train,] %>% dplyr::filter(complete.cases(.))
head(data_train)

# ========================================================
# hSDM model
# ========================================================

# Spatial cells for spatial-autocorrelation
neighborhood <- far$cellneigh_ctry(raster="data/model/fordefor2010.tif",
                                   vector="data/mada/mada38s.shp",
                                   csize=10L, rank=1L)
nneigh <- neighborhood[[1]]
adj <- neighborhood[[2]]
cell_in <- neighborhood[[3]]
ncell <- neighborhood[[4]]

# Udpate cell number in dataset
cell_rank <- vector()
for (i in 1:nrow(data_train)) {
  cell_rank[i] <- which(cell_in==data_train$cell[i])-1 # ! cells start at zero
}
data_train$cell <- cell_rank

# Formula
data_train$trials <- 1  # Set number of trials to one
formula <- paste0("I(1-fordefor2010) + trials ~ C(sapm) + scale(altitude) + scale(slope) +",
                  "scale(dist_defor) + np.power(scale(dist_defor),2) + ",
                  "scale(dist_edge) + ",
                  "scale(dist_road) + scale(dist_town) + cell")

# Model
mod_binomial_iCAR <- far$model_binomial_iCAR(
  # Observations
  suitability_formula=formula, data=r_to_py(data_train),
  # Spatial structure
  n_neighbors=np_array(nneigh,dtype="int32"), neighbors=np_array(adj,dtype="int32"),
  # Environment
  eval_env=-1L,
  # Chains
  burnin=1000L, mcmc=1000L, thin=1L,
  # Starting values
  beta_start=-99)

# # Save model specifications
# model <- list(_x_design_info=mod_binomial_iCAR$"_x_design_info",
# 							 betas=mod_binomial_iCAR$betas)
# py_save_object(model, "output/model.pkl",
#  							 pickle = "pickle")

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
rho <- rep(-9999,ncell)  # -9999 will be considered as nodata
rho[cell_in+1] <- mod_binomial_iCAR$rho

# Resample
far$resample_rho(rho=r_to_py(rho), input_raster="data/model/fordefor2010.tif",
                 output_file="output/rho.tif",
                 csize_orig=10L, csize_new=1L)

# Plot random effects
far$plot$rho("output/rho_orig.tif",output_file="output/rho_orig.png")
far$plot$rho("output/rho.tif",output_file="output/rho.png")

# Plot with R
mada <- rgdal::readOGR(dsn="data/mada",layer="mada38s")
r.rho_orig <- raster("output/rho_orig.tif")
r.rho <- raster("output/rho.tif")
rho_plot(r.rho_orig, mada, output_file="output/rho_orig_ggplot.png",
         quantiles_legend=c(0.025,0.975),width=4.5, height=8)
rho_plot(r.rho, mada, output_file="output/rho_ggplot.png",
         quantiles_legend=c(0.025,0.975),width=4.5, height=8)

# ========================================================
# Predicting spatial probability of deforestation
# ========================================================

# Compute predictions
far$predict_raster_binomial_iCAR(mod_binomial_iCAR, var_dir="data/model",
                                 input_cell_raster="output/rho.tif",
                                 input_forest_raster="data/model/fordefor2010.tif",
                                 output_file="output/pred_binomial_iCAR.tif",
                                 blk_rows=128L)

# Plot predictions
far$plot$prob("output/pred_binomial_iCAR.tif",output_file="output/pred_binomial_iCAR.png")

# ========================================================
# Predicting forest cover
# ========================================================

deforest <- far$deforest(input_raster="output/pred_binomial_iCAR.tif",
                         hectares=4000000,
                         output_file="output/forest_cover_2050.tif",
                         blk_rows=128L)

# Plot future forest cover
far$plot$fcc("output/forest_cover_2050.tif",output_file="output/forest_cover_2050.png")

# ========================================================
# Model performance using validation data-set
# ========================================================

source("R/perf.R")

# Udpate cell number in dataset
cell_rank <- vector()
for (i in 1:nrow(data_valid)) {
  cell_rank[i] <- which(cell_in==data_valid$cell[i])-1 # ! cells start at zero
}
data_valid$cell <- cell_rank

# Performance iCAR
set.seed(1234)
performance_iCAR <- performance_index(data_valid,mod_binomial_iCAR)
performance_nochange <- performance_index(data_valid,"nochange")
performance_null <- performance_index(data_valid,"null")

# ========================================================
# Model comparison
# ========================================================

# Full model
deviance_full <- 0

# Null model
formula_null <- "I(1-fordefor2010) ~ 1"
dmat_null <- patsy$dmatrices(formula_null, data=r_to_py(data_train), NA_action="drop",
                             return_type="dataframe", eval_env=-1L)
Y <- dmat_null[[1]]
X_null <- dmat_null[[2]]
mod_null <- sm$GLM(Y, X_null, family=sm$families$Binomial())$fit()
deviance_null <- mod_null$deviance

# Model nsre with no spatial random effects
formula_nsre <- paste0("I(1-fordefor2010) ~ C(sapm) + scale(altitude) + ",
                       "scale(slope) + scale(dist_defor) + I(scale(dist_defor)*scale(dist_defor)) + ",
                       "scale(dist_edge) + scale(dist_road) + scale(dist_town)")
mod_nsre <- smf$glm(formula_nsre, r_to_py(data_train), family=sm$families$Binomial(), eval_env=-1L)$fit()
deviance_nsre <- mod_nsre$deviance
# Summary nsre
sink(file="output/summary_mod_binomial_nsre.txt")
print(mod_nsre$summary())
sink()
# Performance nsre
performance_nsre <- performance_index(data_valid,mod_nsre)

# Random Forest
source_python("Python/model_random_forest.py")
formula_rf <- paste0("I(1-fordefor2010) ~ sapm + altitude + ",
                     "slope + dist_defor + dist_edge + dist_road + dist_town")
formula_rf_xy <- paste0("I(1-fordefor2010) ~ sapm + altitude + ",
                        "slope + dist_defor + dist_edge + dist_road + dist_town + X + Y")
mod_rf <- model_random_forest(formula=formula_rf, data=data_train, 
                              eval_env=-1L, n_estimators=500L, n_jobs=2L)
mod_rf_xy <- model_random_forest(formula=formula_rf_xy, data=data_train, 
                                 eval_env=-1L, n_estimators=500L, n_jobs=2L)
# Performance rf
performance_rf <- performance_index(data_valid,mod_rf)
performance_rf_xy <- performance_index(data_valid,mod_rf_xy)

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

# =================
# Model differences
# =================

# Compute predictions with logistic regression
fig_pred_nsre <- far$predict_raster(mod_nsre, var_dir="data/model",
                                    input_forest_raster="data/model/fordefor2010.tif",
                                    output_file="output/pred_nsre.tif",
                                    blk_rows=128L, transform=TRUE)

# Plot predictions
far$plot$prob("output/pred_nsre.tif",output_file="output/pred_nsre.png")

# Future forest cover
forest_cover_nsre <- far$deforest(input_raster="output/pred_nsre.tif",
                                  hectares=4000000,
                                  output_file="output/forest_cover_2050_nsre.tif",
                                  blk_rows=128L)

# Plot future forest cover
far$plot$fcc("output/forest_cover_2050_nsre.tif",output_file="output/forest_cover_2050_nsre.png")

# Raster of differences
far$differences(inputA = "output/forest_cover_2050.tif",inputB = "output/forest_cover_2050_nsre.tif",
                output_file = "output/diff_iCAR_nsre_2050.tif", blk_rows = 1L)
# Plot the raster of differences
far$plot$differences("output/diff_iCAR_nsre_2050.tif",
                     borders = "data/mada/mada38s.shp",
                     output_file = "output/diff_iCAR_nsre_2050.png")

# With ggplot2
r_diff <- raster("output/diff_iCAR_nsre_2050.tif")
mada <- rgdal::readOGR(dsn="data/mada",layer="mada38s")
rect_df <- data.frame(xmin=c(346000,793000), xmax=c(439000,886000),
                      ymin=c(7387000,7815000), ymax=c(7480000,7908000),
                      id=c(1,2))
diff_plot(input_raster = r_diff, input_vector = mada, maxpixels=10e5,
          rect = rect_df, output_file = "output/diff_iCAR_nsre_2050_ggplot.png")
diff_plot(input_raster = r_diff, input_vector = mada, maxpixels=10e5,
          ext = extent(346000, 439000, 7387000, 7480000),
          rect = rect_df[1,], output_file = "output/diff_iCAR_nsre_2050_ggplot_zoom1.png")
diff_plot(input_raster = r_diff, input_vector = mada, maxpixels=10e5,
          ext = extent(793000, 886000, 7815000, 7908000),
          rect = rect_df[2,], output_file = "output/diff_iCAR_nsre_2050_ggplot_zoom2.png")

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
far$plot$fcc("output/fcc_2010_2014_obs.tif",output_file="output/fcc_2010_2014_obs.png")

# Number of deforested pixels
count_d <- far$countpix("output/fcc_2010_2014_obs.tif", value=0L, blk_rows=128L)
count_d$area # 392601 ha

# ==============
# Comp with iCAR

# Predictions iCAR 2010-2014
fcc_iCAR_2014 <- far$deforest(input_raster="output/pred_binomial_iCAR.tif",
                              hectares=count_d$area,
                              output_file="output/fcc_2014_iCAR.tif",
                              blk_rows=128L)

# Correlation between maps at ~10km2
val_iCAR_df <- far$validation_npix(r_pred="output/fcc_2014_iCAR.tif",
                                   r_obs="output/fcc_2010_2014_obs.tif",
                                   output_file="output/npix_iCAR.txt", 
                                   value_f=1,
                                   value_d=0,
                                   square_size=33L)

# Tidy dataset
val_iCAR <- val_iCAR_df %>%
  dplyr::mutate(obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f),
                obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d)) %>%
  dplyr::mutate(box0=0:(nrow(val_iCAR_df)-1)) %>%
  dplyr::filter(!(obs_f==0 & obs_d==0)) # Remove squares with no forest

# Coordinates of the centers of the 10-km boxes
e <- extent(raster("output/fcc_2014_iCAR.tif"))
box_size0 <- 30*33
ncol <- ceiling((xmax(e)-xmin(e))/box_size0)
nrow <- ceiling((ymax(e)-ymin(e))/box_size0)
y_box <- floor(val_iCAR$box0/ncol) # quotient
x_box <- val_iCAR$box0 %% ncol # remainder
y <- ymax(e)-box_size0*(y_box+0.5)
x <- xmin(e)+box_size0*(x_box+0.5)
# Box number from x,y coordinates and box size
coeff <- c(1,2,5,10,15,25,50,75,100,150)
box_size <- box_size0 * coeff # box_size in m
CV <- R2 <- Cor <- vector()
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
for (i in 1:length(box_size)) {
  # Box size
  b <- box_size[i]
  # Number of boxes
  ncol <- ceiling((xmax(e) - xmin(e)) / b)
  nrow <- ceiling((ymax(e) - ymin(e)) / b)
  nbox <- ncol * nrow
  # Box identification
  J <- floor((x - xmin(e)) / b)
  I <- floor((ymax(e) - y) / b)
  box <- I * ncol + J
  # Sum deforested areas by box
  obs_d <- val_iCAR %>% mutate(box=box) %>%
    group_by(box) %>% summarise(sum_obs_d_ha=sum(obs_d_ha)) %>%
    pull(sum_obs_d_ha)
  pred_d <- val_iCAR %>% mutate(box=box) %>%
    group_by(box) %>% summarise(sum_pred_d_ha=sum(pred_d_ha)) %>%
    pull(sum_pred_d_ha)
  #plot(obs_d,pred_d)
  CV[i] <- CV_f(obs_d,pred_d)
  R2[i] <- R2_f(obs_d,pred_d)
  Cor[i] <- cor(obs_d,pred_d,method="pearson")
}
plot(coeff,CV,type="l",xlab="Box size (km)")
plot(coeff,R2,type="l",xlab="Box size (km)")
plot(coeff,Cor,type="l",xlab="Box size (km)")
Cor_iCAR <- Cor

# ==============
# Comp with nsre

# Predictions nsre 2010-2014
fcc_nsre_2014 <- far$deforest(input_raster="output/pred_nsre.tif",
                              hectares=count_d$area,
                              output_file="output/fcc_2014_nsre.tif",
                              blk_rows=128L)

# Correlation between maps at ~10km2
val_nsre_df <- far$validation_npix(r_pred="output/fcc_2014_nsre.tif",
                                   r_obs="output/fcc_2010_2014_obs.tif",
                                   output_file="output/npix_nsre.txt", 
                                   value_f=1,
                                   value_d=0,
                                   square_size=33L)
val_nsre_df <- read.table("output/npix_nsre.txt",header=TRUE,sep=",")

# Tidy dataset
val_nsre <- val_nsre_df %>%
  dplyr::mutate(obs_f_ha=npix2ha(obs_f),pred_f_ha=npix2ha(pred_f),
                obs_d_ha=npix2ha(obs_d),pred_d_ha=npix2ha(pred_d)) %>%
  dplyr::mutate(box0=0:(nrow(val_nsre_df)-1)) %>%
  dplyr::filter(!(obs_f==0 & obs_d==0)) # Remove squares with no forest

# Coordinates of the centers of the 10-km boxes
e <- extent(raster("output/fcc_2014_nsre.tif"))
box_size0 <- 30*33
ncol <- ceiling((xmax(e)-xmin(e))/box_size0)
nrow <- ceiling((ymax(e)-ymin(e))/box_size0)
y_box <- floor(val_nsre$box0/ncol) # quotient
x_box <- val_nsre$box0 %% ncol # remainder
y <- ymax(e)-box_size0*(y_box+0.5)
x <- xmin(e)+box_size0*(x_box+0.5)
# Box number from x,y coordinates and box size
coeff <- c(1,2,5,10,15,25,50,75,100,150)
box_size <- box_size0 * coeff # box_size in m
CV <- R2 <- Cor <- vector()
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
for (i in 1:length(box_size)) {
  # Box size
  b <- box_size[i]
  # Number of boxes
  ncol <- ceiling((xmax(e) - xmin(e)) / b)
  nrow <- ceiling((ymax(e) - ymin(e)) / b)
  nbox <- ncol * nrow
  # Box identification
  J <- floor((x - xmin(e)) / b)
  I <- floor((ymax(e) - y) / b)
  box <- I * ncol + J
  # Sum deforested areas by box
  obs_d <- val_nsre %>% mutate(box=box) %>%
    group_by(box) %>% summarise(sum_obs_d_ha=sum(obs_d_ha)) %>%
    pull(sum_obs_d_ha)
  pred_d <- val_nsre %>% mutate(box=box) %>%
    group_by(box) %>% summarise(sum_pred_d_ha=sum(pred_d_ha)) %>%
    pull(sum_pred_d_ha)
  #plot(obs_d,pred_d)
  CV[i] <- CV_f(obs_d,pred_d)
  R2[i] <- R2_f(obs_d,pred_d)
  Cor[i] <- cor(obs_d,pred_d)
}
plot(coeff,CV,type="l",xlab="Box size (km)")
plot(coeff,R2,type="l",xlab="Box size (km)")
plot(coeff,Cor,type="l",xlab="Box size (km)")
Cor_nsre <- Cor

# Compare correlation
plot(coeff,Cor_iCAR,type="l",xlab="Box size (km)",ylim=c(0,0.35))
lines(coeff,Cor_nsre)

# Backup correlation estimates
corr_df <- data.frame(model=c("nsre","iCAR"),
                      corr=c(corr_nsre,corr_iCAR))
write.table(corr_df,"output/corr_df.txt",row.names=FALSE,sep=",")

# ========================================================
# End
# ========================================================