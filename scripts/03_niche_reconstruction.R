#This script takes the true and jittered sample presence points, calibrates a PCA on the true sample and then compares the centroids of the sampled points to the true models centroid. The primary output of this script are the niche biplot graphics, the scree plots and the boxplots of centroid distance for each level of temporal uncertainty.
cat("\n***SCRIPT 03_niche_reconstruction.R BEGIN***\n")

# Libraries ---------------------------------------------------------------

#The following package install lines have a cavaet for the OSC, so that I do not install the same packages multiple times unnecessarily

#Which packages need to be installed
req_pkgs <- c("dplyr", "tidyr", "stringr",
              "terra", "tidyterra", "sf",
              "vegan", "doSNOW", "foreach",
              "ggplot2", "viridis", "ggdensity", "ggrepel",
              "tictoc", "beepr")

#Which ones are currently missing
miss_pkgs <- req_pkgs[!req_pkgs %in% row.names(installed.packages())]

#OSC package installer
install.packages(miss_pkgs,
                 repos = "https://cran.case.edu/",
                 verbose = FALSE,
                 quiet = TRUE)

#Data manipulation
library(dplyr)
library(tidyr)
library(stringr)

#Geospatial data
library(terra)
library(tidyterra)

#Multivariate analysis
library(vegan)

#Parallelisation
library(doSNOW)
library(foreach)

#Visualisation
library(ggplot2)
library(viridis)
library(ggdensity)
library(ggrepel)
library(patchwork)

#Time keepers
library(tictoc)
library(beepr)

# True PCA calibration ----------------------------------------------------

#OSC printer
cat("***Deriving niche centroids...")

#Reading in the sampled environmental data
sample_env_df <- read.csv("./output/virtual_species/sample_env_df.csv")

#Creating unique vectors to iterate through in the for loop
spp_periods <- paste(rep(unique(sample_env_df$species), 
                         each = length(unique(sample_env_df$period))), 
                     unique(sample_env_df$period), 
                     sep = ":")
errors <- unique(sample_env_df$error_range)
replicates <- unique(sample_env_df$rep)

#Initialising an empty list
env_fit_ls <- list()

#For loop to calculate the true centroids and the position of each sample centroid relative to these
#Sampling rasters and fitting a PCA to the 'true' environmental sample
set.seed(42); tic("Niche Centroids"); for(i in seq_along(spp_periods)){
  
  cat("\n- Beginning", spp_periods[i], "-\n")
  
  #Initialising an empty dataframe at the start of the loop
  if(i == 1){
    centroid_df <- data.frame()
  }  
  
  #Subsetting the environmental layer summaries to rescale the values to unit variance
  time_spp_df <- filter(sample_env_df,
                        species == str_split_i(spp_periods[i], ":", 1) & period == str_split_i(spp_periods[i], ":", 2))
  
  for(t in seq_along(errors)){
    
    cat("> Time Uncertainty:", errors[t], "...")
    
    #Subsetting down the error window
    time_spp_t_df <- time_spp_df %>%
      filter(error_range == errors[t]) 
    
    #Creating a df to calibrate the PCA
    true_pca_df <- time_spp_t_df %>%
      #Pulling out the first replicate as the values for the true environmental variables should be constant among all replicates
      filter(rep == 1) %>% 
      dplyr::select(contains("_true")) %>% 
      #Cutting out the names for this version so that the PCA has consistent names for the environmental variables with the sample version
      rename_with(.fn = ~ str_remove(.x, pattern = "_true"))
    
    #Creating a PCA of the true sample
    spp_t_env_pca <- prcomp(true_pca_df,
                            center = TRUE, 
                            scale = TRUE)
    
    #Saving the true PCA models for evaluation
    saveRDS(spp_t_env_pca,
            paste0("./output/niche_reconstruction/pca/", 
                   time_spp_df[1,]$species, "_", time_spp_df[1,]$period, 
                   "_true.rds"))
    
    #Calculating the means for each PC axis to determine the centroid
    true_centroid_df <- as.data.frame(spp_t_env_pca$x) %>% 
      summarise(PC1 = mean(PC1),
                PC2 = mean(PC2),
                PC3 = mean(PC3),
                PC4 = mean(PC4),
                PC5 = mean(PC5),
                PC6 = mean(PC6)) %>% 
      mutate(time_period = time_spp_t_df[1,]$period,
             species = time_spp_t_df[1,]$species,
             error_range = 0,
             rep = NA)
    
    #Adding on rows as it iterates through
    centroid_df <- bind_rows(centroid_df, 
                             true_centroid_df)
    
    #Setting up the cluster
    centroid_clust <- snow::makeCluster(3)
    doSNOW::registerDoSNOW(centroid_clust)
    
    #Initialising a list to save each treatment combinations replicates
    centroid_ls <- list()
    
    #Running this portion in parallel
    centroid_ls <- foreach(r = seq_along(replicates),
                           .packages = c("dplyr",
                                         "stringr"))%dopar%{
                             
                             #Iterating through each replicate
                             time_spp_t_r_df <- time_spp_t_df %>%
                               filter(rep == replicates[r]) 
                             
                             #Deriving the position of the sample centroids
                             sample_centroid_df <- time_spp_t_r_df %>%
                               dplyr::select(contains("_sample")) %>% 
                               rename_with(.fn = ~ str_remove(.x, pattern = "_sample")) %>% 
                               na.omit() %>% 
                               predict(object = spp_t_env_pca,
                                       newdata = .) %>% 
                               as.data.frame() %>% 
                               summarise(PC1 = mean(PC1),
                                         PC2 = mean(PC2),
                                         PC3 = mean(PC3),
                                         PC4 = mean(PC4),
                                         PC5 = mean(PC5),
                                         PC6 = mean(PC6)) %>% 
                               mutate(time_period = time_spp_t_r_df[1,]$period,
                                      species = time_spp_t_r_df[1,]$species,
                                      error_range = errors[t],
                                      rep = replicates[r])
                             
                             #Returning the sample centroid dataframe to be added to the list
                             sample_centroid_df
                           }
    
    snow::stopCluster(centroid_clust)
    #Adding on the sample centroids
    centroid_df <- bind_rows(centroid_df, 
                             centroid_ls %>% bind_rows())
    #Iteration tracker
    cat(" completed\n")
  }
};toc();beep()
#Took about _____ seconds with full 100 replicates

#Saving the output in case I mess it up
write.csv(centroid_df, 
          "./output/niche_reconstruction/centroid_df.csv",
          row.names = FALSE)

#Re-reading in the dataframe from the drive
centroid_df <- read.csv("./output/niche_reconstruction/centroid_df.csv")

#OSC printer
cat("... niche centroids completed***\n")
cat("\n***SCRIPT 03_niche_reconstruction.R END***\n")


#Loop to quickly extract the variation explained by each PC
for(i in seq_along(spp_periods)){
  
  cat("\n- Beginning", spp_periods[i], "-\n")
  
  #Subsetting the environmental layer summaries to rescale the values to unit variance
  time_spp_df <- filter(sample_env_df,
                        species == str_split_i(spp_periods[i], ":", 1) & period == str_split_i(spp_periods[i], ":", 2))
  
  #Subsetting down the error window
  time_spp_t_df <- time_spp_df %>%
    filter(error_range == errors[1]) 
  
  #Creating a df to calibrate the PCA
  true_pca_df <- time_spp_t_df %>%
    #Pulling out the first replicate as the values for the true environmental variables should be constant among all replicates
    filter(rep == 1) %>% 
    dplyr::select(contains("_true")) %>% 
    #Cutting out the names for this version so that the PCA has consistent names for the environmental variables with the sample version
    rename_with(.fn = ~ str_remove(.x, pattern = "_true"))
  
  #Creating a PCA of the true sample
  spp_t_env_pca <- prcomp(true_pca_df,
                          center = TRUE, 
                          scale = TRUE)
  
  pc_variation <- spp_t_env_pca$sdev / sum(spp_t_env_pca$sdev)
  
  
  if(i  == 1) {
    
    pca_summ_df <- data.frame(
      species = str_split_i(spp_periods[i], ":", 1),
      time_period = str_split_i(spp_periods[i], ":", 2),
      pc1_var = pc_variation[1],
      pc2_var = pc_variation[2],
      pc3_var = pc_variation[3],
      pc4_var = pc_variation[4],
      pc5_var = pc_variation[5],
      pc6_var = pc_variation[6]
    )}else{
      pca_summ_df <- rbind(pca_summ_df,
                           data.frame(
                             species = str_split_i(spp_periods[i], ":", 1),
                             time_period = str_split_i(spp_periods[i], ":", 2),
                             pc1_var = pc_variation[1],
                             pc2_var = pc_variation[2],
                             pc3_var = pc_variation[3],
                             pc4_var = pc_variation[4],
                             pc5_var = pc_variation[5],
                             pc6_var = pc_variation[6]))
    }
}


#Reading in the true PCA objects to produce screeplots from them
true_pca_ls <- list.files("./output/niche_reconstruction/pca/",
                          full.names = TRUE,
                          pattern = ".rds") %>% 
  lapply(readRDS)

#Naming to keep track of objects
names(true_pca_ls) <- list.files("./output/niche_reconstruction/pca/",
                                 full.names = FALSE,
                                 pattern = "rds") %>% 
  str_remove(pattern = "_true.rds")



#Loop to produce screeplots and save the output
for(p in 1:length(true_pca_ls)){
  
  eigenvalue_df <- data.frame(PC = paste0("PC", 1:6),
                              Variance = (true_pca_ls[[p]]$sdev^2) / sum(true_pca_ls[[p]]$sdev^2) * 100)
  
  identifier_df <- names(true_pca_ls)[p] %>% 
    as.data.frame() %>% 
    mutate(species = str_sub(., start = 1, end = 7),
           species_full = case_when(
             species == "bla_car" ~ "V_Blarina carolinensis",
             species == "lep_cal" ~ "V_Lepus californicus",
             species == "mic_pen" ~ "V_Microtus pennsylvanicus",
             species == "neo_alb" ~ "V_Neotoma albigula"),
           
           time_period = str_sub(., start = 9),
           time_period = factor(time_period,
                                levels = c("LGM", "DG", "HOL")),
           period_name = ifelse(time_period == "HOL", "Holocene",
                                ifelse(time_period == "DG", "Deglacial", "Last Glacial\nMaximum")) %>% 
             factor(levels = c("Last Glacial\nMaximum", "Deglacial", "Holocene")))
  
  if(p == 1){
    
    true_pca_df <- bind_cols(identifier_df, eigenvalue_df) 
  }else{
    true_pca_df <- bind_cols(identifier_df, eigenvalue_df) %>% 
      bind_rows(true_pca_df, .)
  }
   
}


#Summarising the eigenvalue results for the first two PCs
true_pca_df %>% 
  filter(PC %in% c("PC1", "PC2")) %>% 
  group_by(species_full, period_name) %>% 
  summarise(var_1_2_exp = sum(Variance))

#Wide variance explained column
true_pca_df %>% 
  select(species_full, period_name, PC, Variance) %>% 
  pivot_wider(id_cols = c("species_full", "period_name"),
              names_from = "PC", values_from = "Variance")

#Cumulative variance to add labels to the screeplot
pca_cum_sum_df <- true_pca_df %>% 
  select(species_full, period_name, PC, Variance) %>% 
  group_by(species_full, period_name) %>% 
  mutate(csum = cumsum(Variance) %>% 
           round(digits = 1))


#Creating a facetted screeplot of the PCA axes
ggplot(true_pca_df,
       aes(x = reorder(PC, -Variance), 
           y = Variance,
           group = 1)) + 
  
  geom_line(color = "steelblue", 
            linewidth = 1) +
  
  geom_point(color = "steelblue", 
             size = 3) +
  
  geom_text(data = pca_cum_sum_df,
            aes(x = reorder(PC, -Variance), 
                y = Variance + 5,
                group = 1,
                label = csum)) +
  
  geom_col(fill = "steelblue", alpha = 0.5) +
  
  labs(x = "Principal Components", 
       y = "% Variance Explained") +
  
  facet_grid(species_full ~ period_name) +
  
  theme_bw() + 
  
  theme(strip.text.y = element_text(face = "italic",
                                    size = 12),
        strip.text.x = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
  
ggsave("./graphs/niche_reconstruction/true_pca_screeplot.png",
       width = 14.6, height = 10.5)


# Niche Plot Data Prep --------------------------------------------------------

#Reloading here for easy access
errors <- unique(sample_env_df$error_range)

#Modifying the dataframe for aesthetic presentation
centroid_mod_df <- centroid_df %>%

  mutate(error_range = as.character(error_range * 100),
         error_range_names = replace(error_range, error_range == "10000", "LQ"),
         error_range_names = factor(error_range_names,
                                    levels = rev(c("200", "400", "600", "800",
                                                   "1000", "1500", "2500",
                                                   "4000", "LQ"))),
         centroid_dist = abs(PC1) + abs(PC2),
         time_period = factor(time_period,
                              levels = c("LGM", "DG", "HOL")),
         species_full = case_when(
           species == "bla_car" ~ "V_Blarina carolinensis",
           species == "lep_cal" ~ "V_Lepus californicus",
           species == "mic_pen" ~ "V_Microtus pennsylvanicus",
           species == "neo_alb" ~ "V_Neotoma albigula")) %>% 
  mutate(period_name = ifelse(time_period == "HOL", "Holocene",
                              ifelse(time_period == "DG", "Deglacial", "Last Glacial\nMaximum")) %>% 
           factor(levels = c("Last Glacial\nMaximum", "Deglacial", "Holocene")))