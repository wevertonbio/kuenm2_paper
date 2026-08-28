################################################################################
#                                                                              #
# KUENM2 - An R package for detailed calibration and construction of ENM #     #
#                                                                              #
# Authors:   Weverton Trindade    (wevertonf1993@gmail.com)                    #
#            Luis F. Arias-Giraldo (lfarias.giraldo@gmail.com)                 #
#            Luis Osorio-Olvera (luismurao@gmail.com)                          #
#            A. Townsend Peterson (town@ku.edu)                                #
#            Marlon E. Cobos       (manubio13@gmail.com)                       #
#                                                                              #
################################################################################


# Install package from Github
# remotes::install_github("marlonecobos/kuenm2")

# Install package from CRAN
# install.packages("kuenm2")

# Website
browseURL("https://marlonecobos.github.io/kuenm2/index.html")

#Load packages
library(kuenm2)
library(terra)
library(RuHere)
library(mapview)
library(profvis)


#Set working directory (OPTIONAL)
getwd() #Working directory
my_dir <- "my_project"
dir.create(my_dir) #Create a directory

# Create folder to save figures
figure_path <- file.path(my_dir, "Figures")
dir.create(figure_path)

#### Import data ####

# Dataframe with occurrences (longitude and latitude)
data(occ_data, package = "kuenm2") #Data example
head(occ_data)
str(occ_data)

# Visualize occurrences in a interactive map
pts <- spatialize(occ_data, long = "x", lat = "y")
mapview(pts)

#Raster with environmental variables
var <- terra::rast(system.file("extdata", "Current_variables.tif",
                               package = "kuenm2")) #Data example
plot(var)

# Remove cell duplicates
# Excluding records that are not exact coordinate duplicates but are located within the same pixel.
occ_unique <- remove_cell_duplicates(data = occ_data, x = "x", y = "y",
                                     raster_layer = var)

#### Prepare data for model calibration ####
# Check number of not na cells
global(var, fun="notNA")

#Prepare data for maxnet model
?prepare_data
sp_swd <- prepare_data(algorithm = "maxnet", #Algorithm: maxnet or glm
                       occ = occ_unique, #Dataframe with occurrences
                       species = "Myrcia hatschbachii", #Name of the species (OPTIONAL)
                       x = "x", y = "y", #Column names with long and lat
                       raster_variables = var, #SpatRaster of variables
                       n_background = 1000, #Number of backgrounds
                       features = c("l", "q", "lq", "lp", "qp", "lqp"), #Feature classes
                       r_multiplier = c(0.1, 0.5, 1, 2, 3)) #Regularization parameters for maxnet
sp_swd
#Calibration data (Sample with data - SWD)
head(sp_swd$calibration_data)
#Candidate models
head(sp_swd$formula_grid)
#Number of candidate models
nrow(sp_swd$formula_grid)

#### Explore variable distribution for occurrence and background points ####
?explore_calibration_hist
calib_data <- explore_calibration_hist(data = sp_swd, #Output of prepare_data
                                       include_m = TRUE, #Include entire area?
                                       raster_variables = var) #A SpatRaster of variables)

#Plot
plot_calibration_hist(calib_data, lines = TRUE,
                      which_lines = c("range", "cl", "mean"))

# Save plot
png(filename = file.path(figure_path, "1. Histogram.png"), res = 600,
    units = "in", height = 5, width = 8)
plot_calibration_hist(calib_data, lines = TRUE,
                      which_lines = c("range", "cl", "mean"))
dev.off()


#Color legend:
## Gray: entire area
## Blue: background points
## Green: presences
## Solid line: ranges (maximum and minimum)
## Dashed line: confidence interval at 95%
## Dotted line: mean

#### Distribution of data for models ####
# Explore spatial distribution
pbg <- explore_partition_geo(data = sp_swd, raster_variables = var[[1]])

# Plot exploration results in geography
terra::plot(pbg)

# Save plot
png(filename = file.path(figure_path, "2. Partition_geo.png"), res = 600,
    units = "in", height = 6, width = 6.8)
terra::plot(pbg)
dev.off()


#### Similarity assessment in partitions ####
# Run extrapolation risk analysis
mop_partition <- explore_partition_extrapolation(data = sp_swd,
                                                 include_test_background = TRUE)
# Simple plot in geographic space
plot_explore_partition(explore_partition = mop_partition, space = "G",
                       type_of_plot = "simple")
# Simple plot in environmental space
plot_explore_partition(explore_partition = mop_partition, space = "E",
                       type_of_plot = "simple",
                       variables = c("bio_7", "bio_15"))

# Save plots
png(filename = file.path(figure_path, "3.1. Partition_extrapolation_G.png"), res = 600,
    units = "in", height = 1.75, width = 4, pointsize = 6)
plot_explore_partition(explore_partition = mop_partition, space = "G",
                       type_of_plot = "simple")
dev.off()

png(filename = file.path(figure_path, "3.2. Partition_extrapolation_E.png"), res = 600,
    units = "in", height = 1.75, width = 4, pointsize = 6)
plot_explore_partition(explore_partition = mop_partition, space = "E",
                       type_of_plot = "simple",
                       variables = c("bio_7", "bio_15"))
dev.off()



#### Calibrate maxnet models ####
?calibration

# Testing for concave curves
browseURL("https://marlonecobos.github.io/kuenm2/articles/model_calibration.html#concave-curves")

m <- calibration(data = sp_swd, # prepare_data output
                 error_considered = c(10, 20), # values (in %) used to calculate the omission rate.
                 omission_rate = 10, #Omission rate used to select the best models
                 parallel = TRUE, #Run candidate models in parallel?
                 ncores = 8, #Number of cores for parallel processing
                 remove_concave = TRUE, #Whether to test and remove concave curves
                 progress_bar = TRUE) #Whether to show a progress bar
m

# Check time processing and memory usage
time_calibration <- profvis(calibration(data = sp_swd, # prepare_data output
                                        error_considered = c(10, 20), # values (in %) used to calculate the omission rate.
                                        omission_rate = 10, #Omission rate used to select the best models
                                        parallel = TRUE, #Run candidate models in parallel?
                                        ncores = 8, #Number of cores for parallel processing
                                        remove_concave = TRUE, #Whether to test and remove concave curves
                                        progress_bar = TRUE)) #Whether to show a progress bar)
time_calibration


# #See results for all candidate models
# head(m$calibration_results$Summary)
#
# #### Select models that perform the best among all candidates ####
# ?select_models
# # Use another omission rate available in the calibration_results and higher AIC
# colnames(m$calibration_results$Summary)
# m$omission_rate
# m_20 <- select_models(calibration_results = m,
#                      compute_proc = TRUE,
#                      omission_rate = 20)  # New omission rate
#
# # Compare selected models
# m$selected_models
# m_20$selected_models

#### Training partition effects ####
# What happens when a particular partition is removed?
# ID of models that were selected
m$selected_models$ID

# Response curves for model 313
partition_response_curves(calibration_results = m, modelID = 313)

# Save calibration data
saveRDS(m, file.path(my_dir, "calibration_data.rds"))
# Load calibration data
m <- readRDS(file.path(my_dir, "calibration_data.rds"))

#### Fit selected models ####
?kuenm2::fit_selected
m
# Fit
fm <- fit_selected(calibration_results = m, #calibration output
                   n_replicates = 4, #Number of replicates
                   replicate_method = "kfolds") #Replicate type
fm

# See names of selected models
names(fm$Models)

# See models inside Model
names(fm$Models$Model_313)

# Get omission error used to select models and calculate the theshold values
fm$omission_rate

# See Thresholds
fm$thresholds

#### Predict selected models for a single scenario ####
?kuenm2::predict_selected
p <- predict_selected(models = fm, #output of fit_selected
                      new_variables = var) # SpatRaster of variables for predicting


#Plot replicates
# Plot replicates
plot(p[[1]]$Replicates)

# Plot consensus between replicates
plot(p[[1]]$Model_consensus)

# Plot consensus between models
plot(p$General_consensus)

#### Variable response curves ####
# All response curves
all_response_curves(fm)

# Save plot of response curves
png(filename = file.path(figure_path, "5. All_response_curves.png"), res = 600,
    units = "in", height = 3, width = 3, pointsize = 6)
all_response_curves(fm)
dev.off()


# All response curves for model (GAM for variability across replicates)
all_response_curves(fm, modelID = "Model_313", show_variability = TRUE)

# All response curves for model  (add lines for variability across replicates)
all_response_curves(fm, modelID = "Model_313", show_variability = TRUE,
                    show_lines = TRUE)

# Increase area of extrapolation
all_response_curves(fm, extrapolation_factor = 2)

#### Two-Way interaction response plot ####
?kuenm2::bivariate_response

#Response curves (products)
bivariate_response(models = fm, #output of fit_selected
                   modelID = "Model_313", #Which model?
                   variable1 = "bio_1", #Variable in x-axys
                   variable2 = "bio_12") #VAriable in y-axis

bivariate_response(models = fm, #output of fit_selected
                   modelID = "Model_313", #Which model?
                   variable1 = "bio_7", #Variable in x-axys
                   variable2 = "bio_12") #VAriable in y-axis

# Save plot of response curves
png(filename = file.path(figure_path, "6. Bivariate response.png"), res = 600,
    units = "in", height = 3, width = 3.2, pointsize = 6)
bivariate_response(models = fm, #output of fit_selected
                   modelID = "Model_313", #Which model?
                   variable1 = "bio_1", #Variable in x-axys
                   variable2 = "bio_12") #VAriable in y-axis
dev.off()



#### Variable importance ####
?kuenm2::variable_importance
# We can use the var_importance() function to evaluate the relative contribution
# of each variable:
imp_maxnet <- variable_importance(models = fm) #output of fit_selected
# Plot using enmpa package
enmpa::plot_importance(imp_maxnet)

# Save plot of variable importance
png(filename = file.path(figure_path, "7. Variable_importance.png"), res = 600,
    units = "in", height = 2.8, width = 3.2, pointsize = 6)
enmpa::plot_importance(imp_maxnet)
dev.off()

#### Evaluation with independent data ####
# Import independent records (from Neotropic Tree)
data("new_occ", package = "kuenm2")

# See data structure
head(new_occ)

# Check distribution
new_pts <- RuHere::spatialize(new_occ, long = "x", lat = "y")
mapview(pts, col.regions = "forestgreen") +
  mapview(new_pts, col.regions = "firebrick")


# Extract variables to occurrences
new_data <- extract_occurrence_variables(occ = new_occ, x = "x", y = "y",
                                         raster_variables = var)
head(new_data)

# Evaluate models with independent data
res_ind <- independent_evaluation(fitted_models = fm, new_data = new_data)
# Check results
res_ind$evaluation
# All selected models have significant pROC values, but show higher omission rates
# than expected based on the 10% omission threshold used for calibration

# Any records outside the calibration range?
res_ind$predictions$continuous

# Map predictions
new_pts$predictions <- round(res_ind$predictions$continuous$General_consensus.mean, 3)
mapview(new_pts, zcol = "predictions")

# Save the data:
saveRDS(fm,
        file.path(my_dir, "fitted_models.rds"))

# Import data
fm <- readRDS(file.path(my_dir, "fitted_models.rds"))

#### Binarize predictions ####
# Get omission error used to select models and calculate the thesholds
fm$omission_rate

# Get thresholds
thr <- fm$thresholds$consensus
thr

# Binarize models (mean)
plot(p$General_consensus$mean)
mean_bin <- (p$General_consensus$mean >= thr$mean) * 1
plot(mean_bin)
plot(pts, add = TRUE)
plot(new_pts, col = "white", add = TRUE)

#### Detecting changes in predictions ####
# To compare predictions between two single scenarios representing different time periods
# (e.g., present vs. future)

# Read layers representing future conditions
future_var <- terra::rast(system.file("extdata",
                                      "wc2.1_10m_bioc_ACCESS-CM2_ssp126_2081-2100.tif",
                                      package = "kuenm2"))

# Plot future layers
terra::plot(future_var)

# renaming layers to match names of variables used to fit the model
names(future_var) <- sub("bio0", "bio", names(future_var))
names(future_var) <- sub("bio", "bio_", names(future_var))
names(var)
names(future_var)

# Predict under future environmental conditions
p_future <- predict_selected(models = fm,
                             new_variables = future_var,
                             progress_bar = FALSE)

# Plot consensus (mean)
terra::plot(c(p$General_consensus$mean,
              p_future$General_consensus$mean),
            main = c("Present", "Future (SSP 585, 2081-2100)"))

# Compute changes between scenarios
p_changes <- prediction_changes(current_predictions = p$General_consensus$mean, #Spatraster with current predictions
                                new_predictions = p_future$General_consensus$mean, #Spatraster with future (or past) predictions
                                fitted_models = fm, #fitted_model (used to extract thresholds)
                                predicted_to = "future") #Direction of predictions
# Plot result
terra::plot(p_changes)

#### Project Models to Multiple Scenarios ####
# We have a lot of GCMs...
browseURL("https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html")


#### Organize and structure future climate variables from WorldClim ####

# Set the input directory containing the raw future climate variables.
# For this example, the data is located in the "inst/extdata" folder.
in_dir <- system.file("extdata", package = "kuenm2")
list.files(in_dir, pattern = "wc2") #Raw variables downloades from WorldClim

# Create a "_variables" folder in your directory to save results
out_dir_future <- file.path("Future_variables")
dir.create(out_dir_future)

#Use this function to Organize and rename the future climate data, structuring it by year and GCM
?organize_future_worldclim
organize_future_worldclim(input_dir = in_dir, #path to the folder containing the variables
                          output_dir = out_dir_future, #folder where the organized data will be saved
                          overwrite = TRUE, static_variables = var$SoilType)
#Check the folder
list.dirs(out_dir_future)
fs::dir_tree(out_dir_future) #Requires fs package

#### Preparation of data for model projections ####
#We already have the variables of the future organized with organize_future_worldclim
#We also need the raw variables of the present in a folder
# Let's save the present variables (var) in a folder
out_dir_current <- "Current_raw"
dir.create(out_dir_current)
# Save current variables in this directory
writeRaster(var, file.path(out_dir_current, "Variables.tif"), overwrite = T)


#Now, we can prepare projections data using fitted models to check variables
?prepare_projection
pr <- prepare_projection(models = fm, #output of fit_selected
                         present_dir = out_dir_current, #Where the current variable is saved?
                         future_dir = out_dir_future, #Where future variables are saved? (nested)
                         future_period = c("2041-2060", "2081-2100"), #Periods
                         future_pscen = c("ssp126", "ssp585"), #Socioeconomic pathways
                         future_gcm = c("ACCESS-CM2", "MIROC6"), #GCMS
                         past_dir = NULL, #Where past variables are saved?
                         past_period = NULL, #Periods of past (LGM, MID, etc)
                         past_gcm = NULL) #Past GCMS
pr

#Create folder to save projection results
out_dir <- file.path("Projection_results/maxnet")
dir.create(out_dir, recursive = TRUE)

#### Project selected models for multiple scenarios ####
?project_selected
ps <- project_selected(models = fm, #output of fit_selected
                      projection_data = pr, #output of prepare_proj
                      out_dir = out_dir,  #Where to save projections?
                      write_replicates = TRUE,  #write the projections for each replicate?
                      overwrite = TRUE)
#Writing replicates is useful if you want to check the variance coming from differente replicates

#See object
ps
#See files
list.files(out_dir, recursive = TRUE)
fs::dir_tree(out_dir)

#Import some projections
res_present <- rast("Projection_results/maxnet/Present/Present/General_consensus.tif")
res_future <- rast("Projection_results/maxnet/Future/2081-2100/ssp585/MIROC6/General_consensus.tif")
par(mfrow = c(1,2)) #Grid of plot
plot(res_present$mean, main = "Present")
plot(res_future$mean, main = "2081-2100 - SSP85 - MIROC 6")
dev.off()

#### Compute areas of contraction, expansion and stability between scenarios ####
?projection_changes
changes <- projection_changes(model_projections = ps, #output of project_selected
                              write_results = FALSE, #write the raster files containing the computed changes?
                              return_raster = TRUE) #return a list containing all the SpatRasters with the computed changes

# Setting colors
changes_col <- colors_for_changes(changes)

#See results
plot(changes_col$Binarized) #SpatRaster with the binarized models for each GCM
plot(changes_col$Results_by_gcm) #SpatRaster with the computed changes for each GCM
plot(changes_col$Results_by_change$`Future_2041-2060_ssp126`) #List containing the SpatRaster with each computed change for each GCM
plot(changes_col$Summary_changes) #SpatRaster with the general summary

# Save summary changes
png(filename = file.path(figure_path, "8. Summary_changes.png"), res = 600,
    units = "in", height = 3, width = 4, pointsize = 7)
plot(changes_col$Summary_changes) #SpatRaster with the general summary
dev.off()

#### Binarize changes based on the agreement among GCMs ####
?binarize_changes
# Areas of suitability (stable suitable + gain)
future_suitable <- binarize_changes(changes_projections = changes_col,
                                    outcome = "suitable",
                                    n_gcms = 2)
plot(future_suitable)

# Areas of gain
future_gain <- binarize_changes(changes_projections = changes_col,
                                outcome = "gain",
                                n_gcms = 2)
plot(future_gain)

# Areas of loss
future_loss <- binarize_changes(changes_projections = changes_col,
                                outcome = "loss",
                                n_gcms = 2)
plot(future_loss)

# Areas of stability
future_stable <- binarize_changes(changes_projections = changes_col,
                                outcome = "stable-suitable",
                                n_gcms = 2)
plot(future_stable)

# Example with past projections...
browseURL("https://marlonecobos.github.io/kuenm2/articles/projections_chelsa.html")

#### Exploring Model Uncertainty and Variability ####

#### Prediction variance coming from distinct sources in ENMs ####
?projection_variability
v <- projection_variability(model_projections = ps, #output of project_selected
                            write_files = FALSE, #write the raster files containing the computed variance?
                            return_rasters = TRUE) #return a list containing all the SpatRasters with the computed variance

# Variability for the present
terra::plot(v$Present, range = c(0, 0.15))

# Variability for Future_2081-2100_ssp126
terra::plot(v$`Future_2081-2100_ssp126`, range = c(0, 0.25))

# Variability for Future_2081-2100_ssp585
terra::plot(v$`Future_2081-2100_ssp585`, range = c(0, 0.1))

# Save one of the variability results
png(filename = file.path(figure_path, "9. Variability_Future_2081-2100_ssp585.png"), res = 600,
    units = "in", height = 3, width = 3.2, pointsize = 6)
terra::plot(v$`Future_2081-2100_ssp585`, range = c(0, 0.1))
dev.off()

#### Analysis of extrapolation risks using the MOP metric ####
?projection_mop

# maximum annual mean temperature (bio_1) in our model’s calibration data
max(fm$calibration_data$bio_1)

# Import variables from the 2081-2100 period (SSP585, GCM MIROC6)
future_ACCESS_CM2 <- rast(file.path(out_dir_future,
                                    "2081-2100/ssp585/ACCESS-CM2/Variables.tif"))

# Range of values in bio_1
terra::minmax(future_ACCESS_CM2$bio_1)

# Plot
terra::plot(future_ACCESS_CM2$bio_1,
            breaks = c(-Inf, 22.53, Inf))  # Highlight regions with values above 22.53ºC

#Create folder to save MOP results
out_dir_mop <- file.path("MOP_results")
dir.create(out_dir_mop, recursive = TRUE)

# MOP
kmop <- projection_mop(data = fm, projection_data = pr,
                       subset_variables = TRUE,
                       calculate_distance = TRUE,
                       out_dir = out_dir_mop,
                       type = "detailed",
                       overwrite = TRUE, progress_bar = FALSE)
kmop

# Import results from MOP runs
mop_ssp585_2100 <- import_results(projection = kmop,
                                  future_period = "2081-2100",
                                  future_pscen = "ssp585")

# See types of results
names(mop_ssp585_2100)

# Distances
# Large distance values indicate greater dissimilarity from training conditions, highlighting areas with high extrapolation risk
terra::plot(mop_ssp585_2100$distances)

# Basic
# areas where at least one environmental variable is outside training conditions.
terra::plot(mop_ssp585_2100$basic)

# Simple
# number of environmental variables in the projected area that outside training conditions
terra::plot(mop_ssp585_2100$simple)

# Detailed
# Non-analogous conditions towards high values in the ACCESS-CM2 scenario
terra::plot(mop_ssp585_2100$towards_high_end$`Future_2081-2100_ssp585_ACCESS-CM2`)

# Combined


# Sabe plots mops
png(filename = file.path(figure_path, "10. Mop_distance.png"), res = 600,
    units = "in", height = 3, width = 3.5, pointsize = 5.9)
terra::plot(mop_ssp585_2100$distances)
dev.off()

png(filename = file.path(figure_path, "10. Mop_combined_high.png"), res = 600,
    units = "in", height = 3, width = 3.5, pointsize = 5.9)
terra::plot(mop_ssp585_2100$towards_high_combined$`Future_2081-2100_ssp585_ACCESS-CM2`)
dev.off()


# However, when examining the response curves, we find that the only variable showing real extrapolation is bio_7
all_response_curves(fm)
