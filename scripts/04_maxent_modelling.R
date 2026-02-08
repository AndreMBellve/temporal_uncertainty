#This script takes the true and jittered sample presence points, and fits MaxEnt models using dismo on the true and the temporally jittered sample occurrences, and saves the resulting models. Note, the MaxEnt modelling of the jittered samples takes a very long time to run, unless you have access to a high performance computing cluster. Our analysis was originally run on the OSC.
cat("\n***SCRIPT 04_maxent_modelling.R BEGIN***\n")

# Libraries ---------------------------------------------------------------

#The following package install lines have a cavaet for the OSC, so that I do not install the same packages multiple times unnecessarily

#Which packages need to be installed
req_pkgs <- c("dplyr", "stringr", "terra",
              "foreach", "doSNOW",
              "predicts", "tictoc", "beepr")

#Which ones are currently missing
miss_pkgs <- req_pkgs[!req_pkgs %in% row.names(installed.packages())]

#OSC package installer
install.packages(miss_pkgs,
                 repos = "https://cran.case.edu/",
                 verbose = TRUE,
                 quiet = TRUE)

#Data manipulation
library(dplyr)
library(stringr)

#Geospatial data
library(terra)
library(predicts)

#Parallelisation
library(foreach) 
library(doSNOW)

#Time keepers
library(tictoc)
library(beepr)

#rJava specific installation
#Setting the path to the java directory
if(!"rJava" %in% row.names(installed.packages())){

  #Pitzer Java workaround
  Sys.setenv(LDFLAGS = "-L/apps/spack/0.21/pitzer/linux-rhel9-skylake/libiconv/gcc/12.3.0/1.17-bcgrlj2/lib")
  
  install.packages("rJava",
                   repo = "https://cran.case.edu/")
}

#Upping to java memory limit to cope with the larger files when predicting back to the US
library(rJava)

#The following package install lines have a cavaet for the OSC, so that I do not install the same packages multiple times unnecessarily
#Which packages need to be installed
req_pkgs <- c("dismo")

#Which ones are currently missing
miss_pkgs <- req_pkgs[!req_pkgs %in% row.names(installed.packages())]

#OSC package installer
install.packages(miss_pkgs,
                 repos = "https://cran.case.edu/",
                 verbose = TRUE,
                 quiet = TRUE)

#Modelling
library(dismo)

#Upping the memory limit to ensure the models can fit
options(java.parameters = "-Xmx10g")

# Occ & Env data ----------------------------------------------------------

#Reading in the data produced in script 02_sampling_and_uncertainty.r. The spp_env_ls contains the sampled environmental values for presences and background points but not the sampled values. 
spp_env_ls <- readRDS("./output/virtual_species/true_env_pres_abs.rds")

#A nested list of environmental values of presences for each species × sample size replicated r times; Note, these need to be combined with the matching background samples for each period, which are held constant.
sample_ext_ls <- readRDS("./output/virtual_species/sample_env_pres.rds")

# True MaxEnt -------------------------------------------------------------

#Setting up the cluster - 4 spp × 3 time periods
enm_fitting_clust <- snow::makeCluster(12) #There are only 12 models being fitted total
doSNOW::registerDoSNOW(enm_fitting_clust)

#OSC printer
cat(("***MaxEnt modelling loop...\n"))

#This for loop iterates through each of the species × time period, error range and replicate to fit MaxEnt models using dismo. It creates a nested list with the same structure as spp_env_ls with each element being a maxent model. In addition, it stores all of the maxent output files (see ./output/virtual_species/maxent_files/)
set.seed(42); tic("True ENM Fitting")
true_enm_ls <- foreach(i = seq_along(spp_env_ls),
                          .packages = c("dismo",
                                        "dplyr",
                                        "stringr")) %dopar%{

                                          library(rJava)

                                          #Saving the name of the iteration for naming output files
                                          period_spp <- names(spp_env_ls)[i]

                                          #Pulling predictor data
                                          pred_data <- spp_env_ls[[i]] %>%
                                            dplyr::select(contains("_true")) %>%
                                            rename_with(.fn = ~ str_remove(.x, pattern = "_true"))

                                          #Creating presence-dg identifier vector
                                          pres_abs <- spp_env_ls[[i]]$pres_abs

                                          #Creating a directory if need be
                                          model_output <- paste0(getwd(),
                                                                 "/output/virtual_species/maxent_files/",
                                                                 period_spp, "/no_error")

                                          #Checking to see if the model path directory exists
                                          if(!dir.exists(model_output)){
                                            #If not, then it is created
                                            dir.create(model_output,
                                                       recursive = TRUE)
                                          }

                                          #Pulling out the env values for our presences
                                          pres_vals <- spp_env_ls[[i]] %>%
                                            dplyr::select(contains("_true")) %>%
                                            rename_with(.fn = ~ str_remove(.x,
                                                                           pattern = "_true"))

                                          #Fitting the MaxEnt model and saving files to specified path. This model is fitted with five-fold CV. Making a list within the mod_ls list to save the model object (all of this is named)
                                          enm_mod <- dismo::maxent(x = pred_data,
                                                                   p = pres_abs,
                                                                   args = c("-J",
                                                                            "plots=TRUE",
                                                                            "-P"),
                                                                   path = model_output,
                                                                   silent = TRUE)

                                          #Returning the ENM models
                                          enm_mod
                                        }; snow::stopCluster(enm_fitting_clust); toc(); beep()
#Completed all replicates

#Naming indices for clarity and use when writing files
names(true_enm_ls) <- names(spp_env_ls)

#Saving the list containing all maxent models
saveRDS(true_enm_ls,
        "./output/virtual_species/maxent_files/true_maxent_mods.rds")

#OSC printer
cat("completed***\n")

##True ENM Maps -----------------------------------------------------------

#OSC printer
cat("***Creating true distribution maps...")

#Re-reading in the true models to build maps from
true_enm_ls <- readRDS("./output/virtual_species/maxent_files/true_maxent_mods.rds")

#Adding this to see if it fixes the spatraster predictor issue
library(predicts)

#This for loop iterates through each of the species × time period
#Setting up the cluster - 4 spp × 3 time periods
tic("True MaxEnt Dist"); for(i in seq_along(true_enm_ls)){
  # tic("True MaxEnt Dist"); for(i in 5:length(true_enm_ls)){
  #Saving the name of the iteration for naming output files
  period_spp <- names(true_enm_ls)[i]

  #Isolating the period name to read in the correct env rasters for prediction
  period <- str_sub(period_spp, 9)

  #Pulling the correct env for the time period for the prediction
  if(period == "HOL"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-40_.tif",
                           full.names = TRUE) %>%
      rast()
  }
  if(period == "DG"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-115_.tif",
                           full.names = TRUE) %>%
      rast()
  }
  if(period == "LGM"){
    env_pred <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-160_.tif",
                           full.names = TRUE) %>%
      rast()
  }

  #Pulling out just the CV replicate models to iterate through
  true_iter_mods <- true_enm_ls[[i]]

  model_output <- paste0("./output/virtual_species/maxent_files/true_dist/",
                         period_spp, "/")

  #Creating a directory if it doesn't already exist
  if(!dir.exists(model_output)){
    dir.create(model_output,
               recursive = TRUE)
  }

  #Creating heatmap rasters
  predicts::predict(object = true_iter_mods,
                    x = env_pred,
                    filename = paste0(model_output,
                                      period_spp,
                                      ".tif"),
                    overwrite = TRUE)

  #Iteration tracker
  cat("\nCompleted", period_spp)

}; toc(); beep()


#OSC printer
cat("completed***\n")

# Sample MaxEnt -----------------------------------------------------------

#OSC printer
cat("***Sample MaxEnt modelling...\n")

#Initialising the enm storage list using the sample_ext_ls so that it has an identical structure for saving the MaxEnt models
sample_enm_ls <- sample_ext_ls

#This for loop iterates through each of the time periods, data types, species, sample sizes and replicate samples of environmental occurrences and fits a MaxEnt model using dismo. It creates a nested list with the same structure as env_ls with each element being a maxent model. In addition, it saves the individual models as .rds files (see ./output/virtual_species/enms) and stores all of the maxent output files (see ./output/virtual_species/maxent_files)
set.seed(42); tic("Sample ENM fitting"); for(i in seq_along(sample_ext_ls)){
  
  #Tracking the period being modelled
  period_spp <- names(sample_ext_ls)[i]
  cat("\n", "* Beginning", period_spp, "*\n")
  
  #Defensive stop to check that the two lists match
  iter_check_1 <- spp_env_ls[period_spp][[1]] %>% 
    slice_head(n = 1) %>% 
    mutate(id = paste(species, period, sep = "_")) %>% 
    .$id
  #A stop is called if they do not match
  if(period_spp != iter_check_1){
    stop()
  }
  
  #Pulling out the background point env values to be fitted in our models as these will remain constant among all runs
  bg_vals <- spp_env_ls[[period_spp]] %>% 
    filter(pres_abs == 0) %>% 
    dplyr::select(contains("_true")) %>% 
    rename_with(.fn = ~ str_remove(.x, pattern = "_true"))
  
  #Iterating through the degraded and unmodified samples
  for(t in seq_along(sample_ext_ls[[i]])){
    
    #Tracking the iteration
    cat("-", names(sample_ext_ls[[i]][t]), "-")
        
    #Setting up the cluster - 4 spp × 3 time periods
    enm_fitting_clust <- snow::makeCluster(14) 
    doSNOW::registerDoSNOW(enm_fitting_clust) 
    
    #Reinitialising the list to avoid underlying conflicts with rJava being able to acess lockfiles and fluch user preferences
    sample_enm_r_ls <- list()
    
    #foreach loop creating predictions for each species in parallel
    sample_enm_r_ls <- foreach(r = seq_along(sample_ext_ls[[i]][[t]]), 
                               .packages = c("dismo", 
                                             "dplyr",
                                             "stringr")) %dopar%{
                                               
                                               #Completing OSC java setup
                                               # dyn.load("/apps/java/21.0.2/lib/server/libjvm.so")
                                               library(rJava)
                                               
                                               #Defensive stop to check that the two lists match
                                               iter_check_2 <- sample_ext_ls[[period_spp]][[t]][[r]] %>% 
                                                 slice_head(n = 1) %>% 
                                                 mutate(id = paste(species, period, sep = "_")) %>% 
                                                 .$id
                                               #A stop is called if they do not match
                                               if(iter_check_1 != iter_check_2){
                                                 stop()
                                               }
                                               
                                               #Creating a directory if need be
                                               model_output <- paste0("./output/virtual_species/maxent_files/",
                                                                      period_spp, "/",
                                                                      names(sample_ext_ls[[i]][t]),
                                                                      "_error/", 
                                                                      "rep_", r,"/")
                                               
                                               #Checking to see if the model path directory exists
                                               if(!dir.exists(model_output)){
                                                 #If not, then it is created
                                                 dir.create(model_output, 
                                                            recursive = TRUE)
                                               }
                                               
                                               #Pulling out the env values for our presences 
                                               pres_vals <- sample_ext_ls[[i]][[t]][[r]] %>% 
                                                 dplyr::select(contains("_sample")) %>% 
                                                 rename_with(.fn = ~ str_remove(.x, 
                                                                                pattern = "_sample"))
                                               
                                               #Binding the two env variable dfs into one
                                               pred_data <- rbind(pres_vals, bg_vals)
                                               
                                               #Creating presence-dg identifier vector
                                               pres_abs <- c(rep(1, nrow(pres_vals)), 
                                                             rep(0, nrow(bg_vals)))
                                               
                                               #Fitting the MaxEnt model and saving files to specified path. This model is fitted with five-fold CV. Making a list within the mod_ls list to save the model object (all of this is named)
                                               enm_mod <- dismo::maxent(x = pred_data,
                                                                        p = pres_abs,
                                                                        args = c("-J",
                                                                                 "plots=TRUE",
                                                                                 "-P"), 
                                                                        path = model_output)
                                               
                                               #Returning the ENM models
                                               enm_mod
                                             }
        #Completed all replicates
        snow::stopCluster(enm_fitting_clust)
        
        #Creating a path to save the each replicate set to
        model_save_path <- paste0("./output/virtual_species/maxent_models/",
                                  period_spp,"/")
        
        #Creating the directories if they does not exist
        if(!dir.exists(model_save_path)){
          dir.create(model_save_path,
                     recursive = TRUE)
        }
        
        #Saving each replicate set individually to avoid memory issues
        saveRDS(sample_enm_r_ls,
                paste0(model_save_path,
                       names(sample_ext_ls[[i]][t]),
                       "_error.rds"))
        
        #Removing the object in case it creates issues by hanging around and preventing java from flushing
        remove("sample_enm_r_ls")
        
        #Iteration tracker
        cat(" completed! \n")
      }
  #Completed all sample sizes
}; toc(); beep()

#It takes approximately 7 minutes for each model fit when running in parallel (n rep = 3) - all together about 10.5 hours with 3 replicates run in parallel.

#OSC printer
cat("...sample MaxEnt modelling complete\n")
cat("\n***SCRIPT 04_maxent_modelling.R END***\n")
