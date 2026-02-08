#The primary purpose of this script is to create summaries of environmental variability over the temporal uncertainty range, and plot these against the ability of our models to reconstruct niches or distributions. 

# Libraries ---------------------------------------------------------------
#Data manipulation
library(dplyr)
library(dtplyr)
library(stringr)
library(purrr)
library(tidyr)
library(forcats)
library(tibble)
library(janitor)

#Geospatial
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#Graphing
library(patchwork)
library(ggplot2)

#Parallelisation
library(doParallel)
library(parallel)

#Time keeping
library(tictoc)
library(beepr)

# Reading in species raster data ------------------------------------------

#Listing all the relative abudance files for each species × time period
spp_dist_files <- list.files("./output/virtual_species/raw_maps/distribution/",
                             recursive = TRUE,
                             full.names = TRUE,
                             pattern = ".tif")

#Defining how many presences/absences to sample from the layer
n_pres <- 1000

#Initialising a list...
sp_env_sample_ls <- list()

#For loop to iteratively sample the relative abundance layers to create our true occurrences.
set.seed(42); tic("Sampling points"); for(i in seq_along(spp_dist_files)){
  
  #Cleaning up the file name to get the time period
  time_period <- str_sub(spp_dist_files[i], start = -7, end = -5) %>% 
    str_remove("_")
  
  #Reading in the respective raster file
  sp_rast <- rast(spp_dist_files[i])
  
  #Tracking iteration
  cat("\n", i, "- Beginning", names(sp_rast), time_period, "...")
  
  #Sampling occurrences from the raster using the probabilites as weightings
  sp_sample_df <- spatSample(x = sp_rast,
                              size = n_pres,
                              method = "weights",
                              replace = FALSE, 
                              cells = FALSE,
                              xy = TRUE) %>% 
    #Adding on identifier columns
    mutate(species = names(sp_rast),
           period = time_period,
           mid_point = ifelse(period == "HOL", -40,
                              ifelse(period == "DG", -115, 
                                     -160))) %>% 
    
    #Cleaning up column names for consistency across all dataframes
    dplyr::rename("probability" = names(sp_rast))
  
  #Saving the sample dataframe to our list element
  sp_env_sample_ls[[paste(sp_sample_df[1,]$species, 
                          sp_sample_df[1,]$period,
                          sep = "_")]] <- sp_sample_df
  
  #Iteration tracker
  cat("completed.")
  
}; toc()

#Roughly 10 minutes to run
saveRDS(sp_env_sample_ls, 
        "./output/virtual_species/env_variability/spp_pres_env_sample.rds")

#Re-reading to save having to run the above code
sp_env_sample_ls <- readRDS( "./output/virtual_species/env_variability/spp_pres_env_sample.rds")

#Creating a single dataframe for plotting
sp_env_sample_df <- bind_rows(sp_env_sample_ls)

# Get country-level data from Natural Earth
world <- ne_countries(scale = "medium", 
                      returnclass = "sf")

# Filter for North America
north_america <- world[world$continent == "North America", ]

# Plot to check that the points follow the distribution reasonably well
ggplot() +
  
  geom_sf(data = north_america,
          fill = NA, 
          color = "black") +
  
  geom_point(data = sp_env_sample_df,
          aes(y = y,
              x = x,
              colour = mid_point)) +
  
  coord_sf(xlim = c(-170, -50), 
           ylim = c(5, 85), 
           expand = FALSE) +
  
  facet_grid(species ~ mid_point) + 
  
  theme_minimal()

#Sample points look good - all species are following similar distributions with the Holocene ranges more Northerly than the LGM



# Extracting env values ---------------------------------------------------

#Pulling the names of the environmental variables to match them up later
env_variables <- list.files("./data/rasters/env_predictors/",
                            pattern = "_-40_.tif",
                            full.names = TRUE) %>% 
  rast() %>% 
  names()

#Calculating environmental instability around the time period mid-points
mid_points <- c(-40, -115, -160)
time_ranges <- c(2, 4, 6, 8, 10, 15, 25, 40, 100)

#This for loop iterates through our sample points for each time period * species * error range and extracts the jittered env values before calculating summary metrics of how this jittering relates to variance.
tic("Env value extraction");for(i in seq_along(sp_env_sample_ls)){
  
  #Initialising the df to store results if this is the first iteration
  if(i == 1){
    spp_env_summ_df <- data.frame()
  }
  
  #Storing the current species name from the list to append to the df later
  species_name <- str_sub(names(sp_env_sample_ls)[i], 
                          start = 1, 
                          end = 7)
  
  #Assigning the current list element to make it cleaner coding later
  spp_mid_df <- sp_env_sample_ls[[i]]
  
  
  cat("Beginning", species_name, "|", spp_mid_df[1,]$period)
  #Extracting the 'true' values for each of our sample points from the datas midpoint
  env_mid_df <- list.files("./data/rasters/env_predictors/",
                           pattern = paste0(spp_mid_df[1,]$mid_point, 
                                            "_.tif"),
                           full.names = TRUE) %>%
    #Reading in the desired rasters
    rast() %>% 
    #Extracting the sample point values for the midpoints
    terra::extract(x =., 
                   y = spp_mid_df$cell) %>% 
    #Adding on cells as an identifier - the cells arg in extract not working
    cbind(., cell_id = spp_mid_df$cell) %>% 
    #Pivoting to bind these to the df with all the jittered values later
    pivot_longer(cols = -last_col(),
                 names_to = "variable",
                 values_to = "true_value")
  
  
  #Iterating through the error ranges to create vectors of the time points I need to read in
    for(e in seq_along(time_ranges)){
      
      #This is the standard case because 100 is a place holder for the entirity of the late Quaternary
      if(time_ranges[e] != 100){
        #Creating a vector of all possible years to iterate through
        time_interval <- seq(from = spp_mid_df[1,]$mid_point - time_ranges[e],
                             to = spp_mid_df[1,]$mid_point + time_ranges[e],
                             by = 1) %>% 
          as.character()
      }else{
        #This is the case when the jitter values cover the entirity of the LQ
        time_interval <- seq(from = -200,
                             to = 20) %>% 
          as.character()
      }
      
      #Reading in the rasts within the time intervals
      env_interval_rast <- list.files("./data/rasters/env_predictors/",
                                      pattern = ".tif",
                                      full.names = TRUE) %>% 
        #Subsetting the files that are in my vector of time intervals
        map(time_interval, 
            str_subset, 
            string = .) %>% 
        #Map outputs a set of lists
        unlist() %>%
        #Reading in the desired rasters
        rast()

      #Extracting the environmental values from all of the rasters
      env_range_spp_df <- terra::extract(x = env_interval_rast,
                                         y = spp_mid_df$cell) %>% 
        
        #Adding on cells as an identifier - the cells arg in extract not working
        cbind(., cell_id = spp_mid_df$cell) %>% 
        
        
        pivot_longer(cols = -last_col(),
                     names_to = "variable",
                     values_to = "value") %>% 
        
        na.omit() %>% 
        
        lazy_dt() %>% 
        
        #Binding on the mid points so that I can find the difference between the 'true' env values and the jittered ones. This works because cbind vectorising the bind, making the shorter env_mid_df replicate to match the length of the env_range_spp_df.
        dplyr::left_join(env_mid_df, 
                         by = c("cell_id", "variable")) %>% 
        
        mutate(env_error = value - true_value) %>% 
        
        group_by(variable) %>% 
        
        summarise(num_pts = n(),
                  env_mean = mean(env_error),
                  env_sd = sd(env_error),
                  env_min = min(env_error, na.rm = TRUE),
                  env_q0.05 = quantile(env_error, probs = c(0.05)),
                  env_q0.25 = quantile(env_error, probs = c(0.25)),
                  env_q0.5 = quantile(env_error, probs = c(0.5)),
                  env_q0.75 = quantile(env_error, probs = c(0.75)),
                  env_q0.95 = quantile(env_error, probs = c(0.95)),
                  env_max = max(env_error, na.rm = TRUE),
                  env_rmse = sqrt(mean(env_error^2)),
                  env_mae = median(abs(env_error))) %>% 
        
        as.data.frame() %>% 
        
        mutate(error_range = time_ranges[e],
               mid_point = spp_mid_df[1,]$mid_point,
               spp_code = species_name)
        

      
      #In the darkness bind them
      spp_env_summ_df <- rbind(spp_env_summ_df, 
                               env_range_spp_df)
      
      cat("...", time_ranges[e])
    }
  cat(" finished" ,"\n")
  }; toc()
      
  
write.csv(spp_env_summ_df,
          "./output/environmental_variance/spp_uncertainty_env_values.csv",
          row.names = FALSE)

#Creating summary statistics for variance
spp_env_summ_df <- read.csv("./output/environmental_variance/spp_uncertainty_env_values.csv") %>% 
  mutate(period_name = ifelse(mid_point == -40, "Holocene",
                              ifelse(mid_point == -115, "Deglacial",
                                     "Last Glacial\nMaximum")) %>% 
           factor(levels = c("Holocene", 
                             "Deglacial", 
                             "Last Glacial\nMaximum")),
         
         species_full = fct_recode(spp_code,
                                        "V_Blarina carolinensis" = "bla_car",
                                        "V_Neotoma albigula" = "neo_alb",
                                        "V_Lepus californicus" = "lep_cal",
                                        "V_Microtus pennsylvanicus" = "mic_pen"))


# Min-Max Feature Scaling -------------------------------------------------

#Reading in the sampled environmental data
sample_na_df <- read.csv("./output/virtual_species/sample_env_df.csv") %>% 
  filter(pres_abs == 1) %>% 
  dplyr::select(c(species, period, error_range, rep, ends_with("_sample"))) %>% 
  na.omit() %>% 
  group_by(species, period, error_range, rep) %>% 
  summarise(num_na = 100 - n())

#Reading in prediction accuracy data
pred_comp_df <- read.csv("./output/virtual_species/prediction_accuracy/prediction_accuracy_clean.csv") #%>% 


#Reading in the data on centroids for the sample models and calculating the  mean absolute distance from the 'true' models (this is a Euclidean distance metric)
centroid_dist_df <- read.csv("./output/niche_reconstruction/centroid_df.csv") %>% 
  mutate(error_range = error_range * 100) %>% 
  group_by(time_period, species, error_range) %>% 
  summarise(across(starts_with("PC"),
                   ~sum(.x),
                   .names = "mean_{.col}")) %>% 
  ungroup() %>% 
  filter(error_range != 0) %>% 
  rowwise() %>%
  mutate(mean_abs_dist = mean(abs(c_across(starts_with("mean_"))))) %>%
  ungroup() %>% 
  dplyr::select(!starts_with("mean_PC"))

#TO BE UPDATED TO MIN/MAX AFTER RERUN
spp_env_iqr_df <- spp_env_summ_df %>% 
  
  mutate(env_iqr = env_q0.75 - env_q0.25) %>% 
  
  #Grouping for calculating min-max values for the preceeding mutate call
  group_by(variable, species_full) %>% 
  
  mutate(min_val = min(env_iqr),
         max_val = max(env_iqr),
         range = max_val - min_val) %>%
  
  ungroup() %>% 
  
  #Min-Max Feature Scaling
  mutate(across(starts_with("env_"),
                ~ (. - min_val) / (max_val - min_val),
                .names = "scaled_{.col}")) %>% 
  
  #Joining on permutational importance scores
  left_join(true_perm_df,
            by = join_by(spp_code == species,
                         period_name,
                         variable)) %>%

  group_by(error_range, period_name,
           species_full, spp_code) %>%
  
  summarise(across(starts_with("scaled_"),
                   ~sum(.x, na.rm = TRUE),
                   .names = "sum_{.col}"),
            .groups = "keep") %>%
  mutate(error_range = error_range * 100)
  
  
#Binding on prediction error
pred_plot_df <- pred_comp_df %>% 
  
  distinct(species, time_period, 
           mae, rmse,
           species_full, period_name, 
           error_range) %>% 
  
  #Joining on mean absolute error of sample centroids from true values
  left_join(centroid_dist_df,
            by = join_by(species,
                         time_period,
                         error_range)) %>%
  left_join(spp_env_iqr_df, 
            by = join_by(period_name, 
                         species_full, 
                         error_range)) %>% 
 
  mutate(error_range = as.character(error_range),
         error_range_names = replace(error_range, error_range == "10000", "LQ"),
         error_range_names = factor(error_range_names,
                                    levels = c("200", "400", "600", "800",
                                               "1000", "1500", "2500",
                                               "4000", "LQ")),
         period_name = factor(period_name,
                              levels = c("Holocene", 
                                         "Deglacial", 
                                         "Last Glacial\nMaximum"))) %>% 
  group_by(period_name, species_full, error_range_names,
           across(starts_with("sum_scaled_"))) %>%
  summarise(rmse_mean = mean(rmse, na.rm = TRUE),
            centroid_dist = mean(mean_abs_dist))