#This script loads the  probability of occurrence maps of the true MaxEnt models created in script 04_maxent_modelling.R and the jittered sample MaxEnt models to make comparisons between their predictions and assess the sample models prediction error. The primary output of this script are  the boxplots of prediction for each level of temporal uncertainty.
cat("\n***SCRIPT 05_prediction_accuracy.R BEGIN***\n")

# Libraries ---------------------------------------------------------------

#The following package install lines have a cavaet for the OSC, so that I do not install the same packages multiple times unnecessarily

#Which packages need to be installed
req_pkgs <- c("dplyr", "stringr", "janitor", "terra")

#Which ones are currently missing
miss_pkgs <- req_pkgs[!req_pkgs %in% row.names(installed.packages())]

#OSC package installer
install.packages(miss_pkgs,
                 repos = "https://cran.case.edu/",
                 verbose = FALSE,
                 quiet = TRUE)

#Data manipulation
library(dplyr)
library(stringr)
library(janitor)

#Geospatial data
library(terra)

#rJava specific installation
#Setting the path to the java directory
if(!"rJava" %in% row.names(installed.packages())){
  
  #Pitzer Java workaround
  Sys.setenv(LDFLAGS = "-L/apps/spack/0.21/pitzer/linux-rhel9-skylake/libiconv/gcc/12.3.0/1.17-bcgrlj2/lib")
  
  install.packages("rJava",
                   repo = "https://cran.case.edu/")
}

library(rJava)

#OSC package installer
install.packages(c("dismo", "foreach", "doSNOW", 
                   "ggplot2", "tictoc", "beepr"),
                 repos = "https://cran.case.edu/")

#Modelling
library(dismo)

#Parallelisation
library(foreach) 
library(doSNOW)

#Visualisation
library(ggplot2)

#Time keepers
library(tictoc)
library(beepr)


# True cell sampling ------------------------------------------------------

#OSC printer
cat("***Sampling true cells...\n")

#Reading in all the 'true' distribution predictions
true_vs_rast <- list.files("./output/virtual_species/maxent_files/true_dist",
                           recursive = TRUE,
                           pattern = ".tif",
                           full.names = TRUE) %>% 
  rast()

#renaming the raster files for clarity
names(true_vs_rast) <- list.files("./output/virtual_species/maxent_files/true_dist",
                                  pattern = ".tif") %>% 
  str_remove(".tif")

#Reading in true MaxEnt distributions
tic("True grid cell samples"); set.seed(666); for(l in 1:nlyr(true_vs_rast)){

  #Saving the species name from the raster layers name
  species <-  names(true_vs_rast)[[l]] %>%
    str_sub(start = 1L, end = 7L)

  #Taking the time period of the layer from the layers name
  cur_period <-  names(true_vs_rast)[[l]] %>%
    str_sub(start = 9L)


  #Pulling the correct env for the time period for the prediction
  if(cur_period == "HOL"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-40_.tif",
                           full.names = TRUE) %>%
      rast()
  }
  if(cur_period == "DG"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-115_.tif",
                           full.names = TRUE) %>%
      rast()
  }
  if(cur_period == "LGM"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-160_.tif",
                           full.names = TRUE) %>%
      rast()
  }

  #Tracking progress
  cat("Beginning ", species, " | ", cur_period, "...", sep = "")

  #Creating a sample of cells from the current iterations layer
  sample_df <- true_vs_rast[[l]] %>%

    #Converting to a matrix for efficiency and sampling, and stripping out NaN values
    values(na.rm = TRUE,
           dataframe = TRUE) %>%

    #Converting to a dataframe to bind on other identifying features
    #Rounding as the "accuracy" is unnecessary here and adds computational burden. Values are returned in cell order by "values()"
    data.frame(species = species,
               time_period = cur_period,
               true_prob = round(.[,1], digits = 2),
               cell = as.numeric(row.names(.))) %>%

    #Grouping to pull an equal number of cells for each rounded probability from the raster
    group_by(true_prob) %>%

    #Sampling 100 of each of these probability values at random
    slice_sample(n = 100) %>%

    #Getting rid of redundant columns that won't have a match when this is all combine into one dataframe
    dplyr::select(species, time_period, cell, true_prob) %>%

    #Just in case grouping causes issue down the line...
    ungroup()

  #Extracting the environmental predictor values
  env_prob_df <- terra::extract(x = env_pred,
                                y = sample_df$cell) %>%

    #Adding identifiers back on
    bind_cols(sample_df, .)

  #predicate to to determine if I need to create a new dataframe or just add to the one that I have
  if(l == 1){
    true_pred_df <- env_prob_df
  }else{
    true_pred_df <- bind_rows(true_pred_df,
                              env_prob_df)
  }

  #Tracking progress
  cat("completed \n")

}; toc(); beep()

#Takes about  minutes to run  sec

#OSC printer
cat("...sampling true cells completed***\n")


#Saving and re-reading to save rerunning
write.csv(true_pred_df,
          "./output/virtual_species/prediction_accuracy/true_prediction_df.csv",
          row.names = FALSE)

true_pred_df <- read.csv("./output/virtual_species/prediction_accuracy/true_prediction_df.csv")

# Sample predictions ------------------------------------------------------

#OSC printer
cat("***Sample model predictions...\n")

#Reading in sample models
sample_model_files <- list.files("./output/virtual_species/maxent_models",
                                 pattern = ".rds",
                                 full.names = TRUE,
                                 recursive = TRUE)

#Iterating through models individually because of memory limits
tic("Making predictions for comparison"); for(m in seq_along(sample_model_files)){
  
  #Pulling out the time error code for this iteration
  iter_time_error <- sample_model_files[m] %>% 
    str_extract("[A-Za-z0-9]+(_error)")
  
  #Pulling out the spp identifier to match with the true list
  iter_species <- sample_model_files[m] %>% 
    str_remove("./output/virtual_species/maxent_models/") %>% 
    str_remove(paste0("/", iter_time_error, ".rds")) %>% 
    str_sub(start = 1L, end = 7L)
  
  #Pulling out the period identifier to match with the true list
  iter_period <- sample_model_files[m] %>% 
    str_remove("./output/virtual_species/maxent_models/") %>% 
    str_remove(paste0("/", iter_time_error, ".rds")) %>% 
    str_sub(start = 9L)
  
  #Tracking progress
  cat("\n Beginning", iter_species, iter_period, "\n")
  cat("-", iter_time_error, "-\n")
  
  #Reading in this iterations model
  sample_iter_mod <- readRDS(sample_model_files[m])
  
  #Filtering the true pred to the spp and time period to be predicted
  true_iter_df <- true_pred_df %>% 
    filter(species == iter_species) %>% 
    filter(time_period == iter_period) %>% 
    as.data.frame()
  
  #Creating a dataframe of predictor values for the sample true cells
  true_iter_comp_df <- true_iter_df %>% 
    dplyr::select(elevation, precip_meanAnnual, precip_seasonality,
                  slope, temp_meanAnnual, temp_seasonality) 
  
  #Upping to java memory limit to cope with the larger files when predicting back to the US
  options(java.parameters = "-Xmx16g")
  
  #Intialising comparison dataframe to fill and starting timer
  pred_comp_df <- data.frame()
  
  #Loop to iterate through each of the model objects and predict occurrence probabilities to compare to the "true" occ probabilities, and extracting model test AUC values
  for(r in seq_along(sample_iter_mod)){
    
    #Pulling out a model replicate
    cur_mod <- sample_iter_mod[[r]]
    
    #Tracking progress
    cat( r, "...")
    
    cur_mod_pred <- dismo::predict(object = cur_mod,
                                   x = true_iter_comp_df) %>% 
      data.frame(modelled_prob = .)

    #Summarising prediction values
    cur_comp_df <- cur_mod_pred %>%
    
      #Binding on the true data to the prediction results
      bind_cols(true_iter_df) %>% 
      #Adding on identifiers and cleaning prediction probabilities
      mutate(error_range = iter_time_error,
             rep = r,
             modelled_prob = round(modelled_prob, digits = 2)) %>% 
      #Binding on the model evaluation data for later
      bind_cols(t(cur_mod@results) %>% 
                  clean_names() %>% 
                  as_tibble() %>% 
                  dplyr::select(training_auc, 
                                ends_with("permutation_importance")))
    
    #Conditional dataframe assembly
    if(nrow(pred_comp_df) == 0){
      pred_comp_df <- cur_comp_df
    }else{
      pred_comp_df <- bind_rows(pred_comp_df,
                                cur_comp_df)
    }
  }
  
  #Saving the output dataframes individually to avoid memory issues
  write.csv(pred_comp_df,
            paste0("./output/virtual_species/prediction_accuracy/model_predictions/",
                   iter_species, "_", iter_period, "_", iter_time_error,".csv"),
            row.names = FALSE)
  
}; toc(); beep()

#OSC printer
cat("...sample model predictions completed***\n")
cat("\n***SCRIPT 05_prediction_accuracy.R END***\n")

