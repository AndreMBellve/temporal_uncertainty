#Script to sample the 'true' relative abundance layers to create a set of occurrences and absences for each combination of time period and species. Subsequently, we will jitter the time periods for these occurrences and then sample from the respective environmental layers.

#This script has four outputs: 
#1) spp_points_ls ("./output/virtual_species/true_env_pres_abs.rds") which contains sampled occurrences and background points from the 'true' species layers. This dataset is a list of dataframes that have the ID of the sampled grid cell (cell), the longitude and latitude of the cell (x, y), the probability of the sp occurrence in this cell (probability), whether the cell was sampled as a presence of a background point (pres_abs; 1 being presence, 0 being absence), the six letter spp code (species), and the sampled period for these occurrences (period).
#2) spp_env_ls ("./output/virtual_species/true_env_pres_abs.rds"), which has an identical list structure and contains dataframes with all the columns listed above, as well as the extracted environmental values for the corresponding cells which have the associated environmental variables names, concatenated with "_true" to indicate these are the real env values for that time point.
#3) sample_ext_ls ("./output/virtual_species/sample_env_pres.rds"), which is a nested list containing our jittered dataframes for each temporal uncertainty bound and replicate. The nesting is as follows: i) species × period; ii) temporal uncertainty; iii) replicate. *Each dataframe only contains presences*, as we are keeping our background points constant among our maxent models. These dataframes has all the columns listed above, as well as an identifier for the periods true midpoint (mid_point), the added error (time_error), the combination of these two, i.e., the jittered time point (jitter_point), the error range it was sampled from (error_range), the replicate identifier of the jittering sample (rep), and the environmental values of these jittered points which have the same name as those of the environmental variables, but with "_sample" concatenated onto the ends of the names.
#4) sample_env_df ("./output/virtual_species/sample_env_df.csv") - Unnested version of the above list for simple cross-referencing and extraction. This contains 24 columns. The 'cell' denotes the raster cell id that the sample/env data was drawn from. 'x' and 'y' are the longitude and latitude, respectively, stored in WGS 84 for these cells. The 'probability' column denotes the probability of suitable habitat (0 - 1) for the virtual species in question. 'pres_abs' is a 0/1 binary indicating whether the row corresponds to a presence or absence (read background point), but this data frame only contains presences. 'period' shows the 2/3 letter time period name abbreviation that the sample points were originally drawn from, while 'mid_point' shows the century that the data were drawn from (CE/BCE, corresponding to positive/negative values respectively, with year 0 being approximately 2000 years ago). Note, these values correspond to centuries (e.g., 2 = 200 years, -115 = 11,500 BCE = 13,500 years before present), to match the typology of the CHELSA-TRacE 21k data. 'elevation_true' (metres above mean sea level), 'precip_mean_true' (kg/m2/year), 'precip_seasonality_true' (coefficient of variation), 'slope_true' (degrees), 'temp_meanAnnual_true' (degrees Celcius) and 'temp_seasonality_true' (standard deviation of monthly means), all give the environmental values for the cell from the corresponding time period that the species samples were originally drawn from. 'time_error' indicates the temporal jittering that was added to the observation, although the values correspond to centuries as described above for the 'mid_point'. 'jitter_point' shows the actual layer that the temporally jittered value was drawn from (using CHELSA-TraCE 21k typology). 'error_range' indicates the bounds of temporal uncertainty that could be added to the observation. 'rep' designates the ID of the temporal jittering applied to the sample observations. 'elevation_sample', 'precip_mean_sample', 'precip_seasonality_sample', 'slope_sample', 'temp_meanAnnual_sample' and 'temp_seasonality_sample' show the values drawn from the temporally jittered sample points from their new environmental raster set (e.g., the rasters from the century listed in 'jitter_point').

cat("\n***SCRIPT 02_sampling_and_uncertainty.R BEGIN***\n")

# Libraries ---------------------------------------------------------------

#The following package install lines have a cavaet for the OSC, so that I do not install the same packages multiple times unnecessarily

#Which packages need to be installed
req_pkgs <- c("dplyr", "stringr",
               "terra", "tidyterra", "sf",
               "ggplot2", "tictoc", "beepr")

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

#Geospatial packages
library(terra)
library(tidyterra)
library(sf)

#Visualisation
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)

#Time keeping
library(tictoc)
library(beepr)

# True Samples ------------------------------------------------------------

#OSC printer
cat("***Creating true samples...\n")

#Listing all the relative abudance files for each species × time period
spp_dist_files <- list.files("./output/virtual_species/raw_maps/distribution/",
                             recursive = TRUE,
                             full.names = TRUE,
                             pattern = ".tif")

#Defining how many presences/absences to sample from the layer
n_pres <- 100
n_abs <- 10000

#Initialising a list...
spp_points_ls <- list()

#For loop to iteratively sample the relative abundance layers to create our true occurrences.
set.seed(42); tic("Sampling points"); for(i in seq_along(spp_dist_files)){
  
  #Cleaning up the file name to get the time period
  time_period <- str_sub(spp_dist_files[i], start = -7, end = -5) %>% 
    str_remove("_")
  
  #Reading in the respective raster file
  spp_rast <- rast(spp_dist_files[i])
  
  #Tracking iteration
  cat("\n", i, "- Beginning", names(spp_rast), time_period, "...")
  
  #Sampling occurrences from the raster using the probabilites as weightings
  spp_sample_df <- spatSample(x = spp_rast,
                              size = n_pres,
                              method = "weights",
                              replace = FALSE, 
                              cells = TRUE,
                              xy = TRUE) %>% 
    
    #Sampling absences
    rbind(spatSample(x = spp_rast,
                     size = n_abs,
                     method = "random",
                     replace = FALSE, 
                     na.rm = TRUE,
                     exhaustive = TRUE,
                     cells = TRUE,
                     xy = TRUE)) %>% 
    
    #Adding on identifier columns
    mutate(species = names(spp_rast),
           pres_abs = c(rep(1, n_pres), 
                        rep(0, n_abs)),
           period = time_period,
           mid_point = ifelse(period == "HOL", -40,
                              ifelse(period == "DG", -115, 
                                     -160))) %>% 
    
    #Cleaning up column names for consistency across all dataframes
    dplyr::rename("probability" = names(spp_rast))
  
  #Saving the sample dataframe to our list element
  spp_points_ls[[paste(spp_sample_df[1,]$species, 
                       spp_sample_df[1,]$period,
                       sep = "_")]] <- spp_sample_df
  
  #Iteration tracker
  cat("completed.")
  
}; toc()
#Roughly 50 minutes to run
saveRDS(spp_points_ls, "./output/virtual_species/pres_abs_points.rds")
#Re-reading to save having to run the above code
spp_points_ls <- readRDS("./output/virtual_species/pres_abs_points.rds")

#OSC printer
cat("completed***\n")

# Environmental Values ----------------------------------------------------

#OSC printer
cat("***Extracting env values for true...\n")

#I am intentionally breaking apart what could be a single for loop for clarity and to help isolate potential errors
tic("Extracting env values")
for(i in seq_along(spp_points_ls)){
  
  #Pulling out the iterations dataframe
  spp_sample_df <- spp_points_ls[[i]]
  
  #Tracking iteration
  cat("\n", i, "- Beginning", names(spp_points_ls)[i], "...")
  
  #Time period identifier
  time_period <- spp_sample_df[1,]$period
  
  #Pulling the correct environmental layer
  if(time_period == "HOL"){
    #Holocene
    env_rast <- list.files("./data/rasters/env_predictors/",
                               pattern = "_-40_.tif",
                               full.names = TRUE) %>%
      rast()
  }
  
  if(time_period == "DG"){
    #Deglacial
    env_rast <- list.files("./data/rasters/env_predictors/",
                              pattern = "_-115_.tif",
                              full.names = TRUE) %>%
      rast()
  }

  if(time_period == "LGM"){
    #Last glacial maximum
    env_rast <- list.files("./data/rasters/env_predictors/",
                               pattern = "_-160_.tif",
                               full.names = TRUE) %>%
      rast()
  }
  
  #Extracting environmental values for the species sample
  spp_env_df <- terra::extract(x = env_rast, 
                               y = spp_sample_df$cell) 
    
  
  #Adding on an identifier for the env value names to mark them as the true environmental values
  names(spp_env_df) <- paste0(names(spp_env_df), "_true")
  
  #Adding the env values onto our original sample dataset
  spp_df <- cbind(spp_sample_df, spp_env_df)
  
  #Initialising an empty list if necessary, with the same structure/names as the points list
  if(i == 1){
    spp_env_ls <- list()
    spp_env_ls[[names(spp_points_ls)[1]]] <- spp_points_ls
  }
  
  #Saving the sample dataframe to our list element
  spp_env_ls[[names(spp_points_ls)[i]]]<- spp_df
  
  #Iteration tracker
  cat("completed.")
};toc();beep()
#Run time is about a minute (3 reps)
saveRDS(spp_env_ls,
        "./output/virtual_species/true_env_pres_abs.rds")
#Re-reading to save having to run the above code
spp_env_ls <- readRDS("./output/virtual_species/true_env_pres_abs.rds")

#OSC printer
cat("completed***\n")

# Uncertainty Jittering ----------------------------------------------------

#OSC printer
cat("***Uncertainty Jittering...")

#Temporal resolution windows
#These values are from our notes - the NA indicates when we sample randomly from the last 23,000 years
time_ranges <- c(2, 4, 6, 8, 10, 15, 25, 40, 100)

#Number of replicates of each level of uncertainty
n_reps <- 100

#Initialising the sample list
sample_env_ls <- list()

#For loop to add temporal uncertainty
set.seed(42); tic("Jitter Loop"); for(i in seq_along(spp_env_ls)){
  
  #Subsetting the true species list
  true_spp_df <- spp_env_ls[[i]] %>% 
    #Tidying up the values to save on space and process more efficiently
    mutate(across(ends_with("_true") | ends_with("probability"), ~ round(.x, digits = 2)),
           mid_point = ifelse(period == "HOL", -40,
                              ifelse(period == "DG", -115, 
                                     -160))) %>% 
    #Filtering down to just the presences, as the absence values are going to be held constant - THIS ASSUMPTION TO BE CHECKED!
    filter(pres_abs == 1)
  
  #Trimming row names because eeewwwwwwww.....
  row.names(true_spp_df) <- NULL
  
  #Initialising a list in a list!
  sample_env_ls[[names(spp_env_ls)[i]]] <- list()
  
  #Cycling through our time ranges of uncertainty
  for(t in seq_along(time_ranges)){
    
    if(time_ranges[t] != 100){
      
      #Creating a sequence of all possible numbers between in our error range to uniformly sample from.
      error_range <- seq(from = -1 * time_ranges[t],
                         to = time_ranges[t],
                         by = 1)
      
      for(r in 1:n_reps){
        
        sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges[t])]][[r]] <- true_spp_df %>% 
          mutate(time_error =  sample(error_range,
                                      size = n_pres,
                                      replace = TRUE),
                 jitter_point = mid_point + time_error,
                 error_range = time_ranges[t],
                 rep = r)
      }
      
      #Exception needed for late Quaternary sampling as it is not possible to simply add error bounds since they are assymetrical relative to the midpoints. Instead of sampling errors, we sample years directly and then back calculate the error 
    }else{
      
      #These are the years being sampled instead of errors to be added around the midpoints like above. This is because the intervals around the points cannot be completely even
      error_range <- seq(from = -200,
                         to = 20)
      
      for(r in 1:n_reps){
        
        sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges[t])]][[r]] <- true_spp_df %>% 
          mutate(time_error = sample(error_range,
                                     size = n_pres,
                                     replace = TRUE),
                 jitter_point = time_error,
                 error_range = time_ranges[t],
                 rep = r) %>% 
          #Rewriting time_error so it is in the same format as the values in our other dfs
          mutate(time_error = jitter_point - mid_point)
      }
    }
  }
}; toc(); beep()
#Takes seconds...(3 reps)

#Saving the output
saveRDS(sample_env_ls,
        "./output/virtual_species/spp_pres_jitter.rds")

#OSC printer
cat("completed***\n")


# Uncertainty extraction --------------------------------------------------

#OSC printer
cat("***Uncertainty sampling...\n")

#Creating one giant dataframe that I will filter through and extract environmental values for. This should be a slightly more efficient way of pulling out the env values as I will not need to load and extract from the same raster multiple times
sample_env_df <- do.call(rbind, 
                         unlist(unlist(sample_env_ls, 
                                       recursive = F), 
                                recursive = F))

#Creating a vector of unique years to iteratively read in and extract values for
year_samples <- unique(sample_env_df$jitter_point)

#For loop to iteratively pull the necessary jittered layers and then extract the env values from them
tic("Extracting env val"); for(i in seq_along(year_samples)){
  
  #Reading in all layers for a given year
  env_rast <- list.files("./data/rasters/env_predictors/",
                         pattern = paste0("_", year_samples[i], "_.tif"),
                         full.names = TRUE) %>%
    rast()
  
  #Filtering the full dataset down to the years currently being iterated through
  year_filter_df <- filter(sample_env_df,
                           jitter_point == year_samples[i])
  
  #Extracting the values from these years - NOTE: This step introduces some NAs as the glaciers move and may cover locations for any particular time period. This is most likely to happen to our Holcene samples, but it appears to be relatively small numbers of points.
  year_extract_df <- terra::extract(x = env_rast, 
                                    y = year_filter_df$cell)
  
  #Modifying the variable name so that I can retain the "true values for the occurrence in the same row
  names(year_extract_df) <- paste0(names(year_extract_df), "_sample")
  
  #Binding the extracted values with the filtered dataset so it has the identifiers we need
  year_extract_df <- bind_cols(dplyr::select(year_filter_df, 
                                             c(cell, x, y, 
                                               species, 
                                               jitter_point, mid_point, 
                                               error_range, rep)),
                               year_extract_df)
  
  #Initialising or adding to a large df that will contain everything to later be separated. This could create problems if it gets too big!!
  if(i == 1){
    sample_extract_df <- year_extract_df
  }else{
    sample_extract_df <- rbind(sample_extract_df, 
                               year_extract_df)
  }
  
}; toc(); beep()
#Took approximately 50 seconds to run (3 reps)

#Joining on the cells so that the extracted values match up with the cells they were extracted from. Including all the identifiers for an added bit of certainty that nothing gets double matched incorrectly!
sample_env_df <- left_join(sample_env_df,
                           sample_extract_df,
                           by = join_by(cell, x, y, 
                                        species, 
                                        jitter_point, mid_point, 
                                        error_range, rep))

#Saving the output for the 03_niche_reconstruction.R script
write.csv(sample_env_df,
          "./output/virtual_species/sample_env_df.csv",
          row.names = FALSE)

#Quick check to see if any values are missing - If TRUE, check to see why these points are missing (likely glaciers & lakes)
any(is.na(sample_env_df$elevation_sample))

#Initialising list with the same structure
sample_ext_ls <- sample_env_ls

#Ugly loop to split it back into a nested list
tic("Sorting df to ls"); for(i in seq_along(sample_env_ls)){
  for(t in seq_along(sample_env_ls[[i]])){
      for(r in seq_along(sample_env_ls[[i]][[t]])){
        sample_ext_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges)[t]]][[r]] <- sample_env_df %>% 
          filter(species == sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges)[t]]][[r]][1,]$species) %>% 
          filter(period == sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges)[t]]][[r]][1,]$period) %>% 
          filter(error_range == sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges)[t]]][[r]][1,]$error_range) %>% 
          filter(rep == sample_env_ls[[names(spp_env_ls)[i]]][[as.character(time_ranges)[t]]][[r]][1,]$rep)
    }
  }  
}; toc(); beep()
#Approximately 10 seconds

#Saving a copy of the sampled presence points
saveRDS(sample_ext_ls,
        "./output/virtual_species/sample_env_pres.rds")

#OSC printer
cat("completed***\n")
cat("\n***SCRIPT 02_sampling_and_uncertainty.R END***\n")

