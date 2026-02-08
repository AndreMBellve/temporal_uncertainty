#Script to load the virtual species and create 'true' distributions for each time period and species

# Libraries ---------------------------------------------------------------

#Data manipulation
library(dplyr)
library(stringr)

#Geospatial packages
library(terra)
library(tidyterra)

#Virtual species creation
library(virtualspecies)

#Mathematical functions
library(EnvStats) 
library(stats)

#Visualisation
library(ggplot2)

#Time keeping
library(tictoc)
library(beepr)

# Functions ---------------------------------------------------------------

#An updated version of the virtualspecies::generateSpFromFun() function adapted to use terra instead of the raster package - this significantly reduces computational time, although further streamlining is definitely possible!
source("./scripts/functions/gen_sp_fun.R")

#Script that loads custom mathematical functions for creating virtual species response curves
source("./scripts/functions/sigmoid_functions.R")

#Creating a list of the response functions to feed into the virtualspecies functions to generate the response curves at each time point - this is a subset of the species we used for the virtual taphonomy study.
source("./scripts/functions/vs_resp_ls.R")

# Paleo-env data ----------------------------------------------------------

#Reading in the environmental data for our three 'true' time periods

#Holocene
hol_mid_rast <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-40_.tif",
                           full.names = TRUE) %>%
  rast()

#Deglacial
dg_mid_rast <- list.files("./data/rasters/env_predictors/",
                          pattern = "_-115_.tif",
                          full.names = TRUE) %>%
  rast()

#Last glacial maximum
lgm_mid_rast <- list.files("./data/rasters/env_predictors/",
                           pattern = "_-160_.tif",
                           full.names = TRUE) %>%
  rast()


# Paleo-distributions -----------------------------------------------------

#List of species names to iterate through
vs_spp <- unique(names(vs_resp_ls))
#Periods to cover
periods <- c("HOL", "DG", "LGM")

#Creating a list to fill with the output of the gen_sp_fun function for each species × time period
vs_paleo_ls <- list()

#For loop to iterate through the species × time period
for(sp in seq_along(vs_spp)){

  #Keeping track of the iteration being progressed
  cat(paste0("Creating ", vs_spp[sp], " modelling... \n"))

  #Creating a list within the list for each species
  vs_paleo_ls[[vs_spp[sp]]] <- list()

  for(p in seq_along(periods)){

    cat(paste0("Beginning ", periods[p], "...\n"))

    #Pulling the correct paleo data for the current period
    if(periods[p] == "HOL"){
      paleo_env_rast <- hol_mid_rast
    }
    if(periods[p] == "DG"){
      paleo_env_rast <- dg_mid_rast
    }
    if(periods[p] == "LGM"){
      paleo_env_rast <- lgm_mid_rast
    }

    #Creating the habitat suitability raster for the current species × period iteration
    vs_paleo_ls[[vs_spp[sp]]][[periods[p]]] <- gen_sp_fun(env.rast = paleo_env_rast,
                                                          parameters = vs_resp_ls[[vs_spp[sp]]],
                                                          rescale = TRUE,
                                                          formula = NULL,
                                                          species.type = "multiplicative",
                                                          rescale.each.response = TRUE,
                                                          plot = FALSE)

    cat(" Finished \n")
  }

}

#Saving the species × period rasters
for(sp in seq_along(vs_spp)){

  for(p in seq_along(periods)){

    #Naming the raster file with just the species name
    names(vs_paleo_ls[[vs_spp[sp]]][[periods[p]]]$suitab.raster) <- vs_spp[sp]

    #Naming the tif file with the sp × period combination
    file_path <- paste0("./output/virtual_species/raw_maps/distribution/", vs_spp[sp], "/",
                        vs_spp[sp], "_", periods[p], ".tif")

    #Writing the raster straight to disc
    writeRaster(vs_paleo_ls[[vs_spp[sp]]][[periods[p]]]$suitab.raster,
                filename = file_path,
                overwrite = TRUE,
                gdal = c("COMPRESS=DEFLATE"))

  }
}

#Saving the vs_paleo_ls output as these contain the rescaled response variables
saveRDS(vs_paleo_ls,
        "./output/virtual_species/raw_maps/true_vs_niches.rds")
#Read in the saved curve
vs_paleo_ls <- readRDS("./output/virtual_species/raw_maps/true_vs_niches.rds")
