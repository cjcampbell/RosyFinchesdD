source('~/PESU_migration_project_directory/PESU_migration/R/00_Setup.R')
source('~/PESU_migration_project_directory/PESU_migration/R/01_loadIsotopeData.R')

# Load direction stats.
if(!exists("distDir_ORSim") ) distDir_ORSim <- readRDS( file.path(wd$bin, "distDir_ORSim.rds") )
dirSummary <- distDir_ORSim %>% 
  dplyr::select(-Sex) %>% 
  dplyr::mutate(theta_rounded = round(theta_from_origin, 0)) %>% 
  group_by(ID, theta_rounded, .drop = F) %>% 
  dplyr::summarise(
    mean_rawVal = mean(value, na.rm = T),
    mean_OR = mean(OR,  na.rm = T )
    ) %>% 
  ungroup()

dir_deets <- dirSummary %>% 
  dplyr::select(-mean_OR) %>% 
  pivot_wider(
    names_from = theta_rounded,
    names_prefix = "probOfOriginAtAngle_" ,
    values_from = mean_rawVal
    )

dir_deets_SHORTER <- distDir_ORSim %>% 
  dplyr::mutate(
    direction = case_when(
      theta_from_origin >= -45 & theta_from_origin < 45 ~ "North",
      theta_from_origin >= 45 & theta_from_origin < 135 ~ "East",
      theta_from_origin >= -135 & theta_from_origin < -45 ~ "West",
      TRUE ~ "South"
    )
    ) %>% 
  dplyr::group_by(ID, direction) %>% 
  dplyr::summarise(
    ave_prob = mean(value, na.rm = T)
    ) %>% 
  pivot_wider(
    values_from = ave_prob,
    names_from = direction
  )

# Load distance stats.
if(!exists("minDistDeets") ) minDistDeets <- readRDS( file.path(wd$bin, "minDistDeets.rds") )
dist_deets <- minDistDeets %>% 
  dplyr::select(-Sex) %>% 
  dplyr::select(ID, threshold, minDist) %>% 
  dplyr::filter(threshold %in% c( 0.25, 0.32, 0.5, 0.68, 0.75) ) %>% 
  arrange(ID, threshold) %>% 
  pivot_wider(
    names_from = threshold,
    names_prefix = "probDistTraveledAtThreshold_" ,
    values_from = minDist
  )

# wns results
if(!exists("subsampledOrigins_withMinDistance") ) subsampledOrigins_withMinDistance <- readRDS( file.path(wd$bin, "subsampledOrigins_withMinDistance.rds") )
wns_deets <- subsampledOrigins_withMinDistance %>% 
  dplyr::select(-Sex) %>% 
  dplyr::select(ID, distFromLineToNearestPdPositiveCounty) %>% 
  dplyr::group_by(ID) %>% 
  group_by(ID) %>% 
  dplyr::summarise(
    
    propSimulationsPassingWithin5kmOfPdPositiveCounty = round( sum(distFromLineToNearestPdPositiveCounty <= 5) / n() , 2 ),
    propSimulationsPassingWithin25kmOfPdPositiveCounty = round( sum(distFromLineToNearestPdPositiveCounty <= 25) / n() , 2 ),
    propSimulationsPassingWithin100kmOfPdPositiveCounty = round( sum(distFromLineToNearestPdPositiveCounty <= 100) / n() , 2 )
    
  ) 


# Combine -----------------------------------------------------------------
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

PESUResults <- mydata_transformed %>% 
  dplyr::select(
    ID, verbatimID, contains("decimal"), contains("meters"),
    everything()
    ) %>% 
  full_join(dist_deets) %>% 
  full_join(wns_deets) %>% 
  full_join(dir_deets_SHORTER) %>% 
  full_join(dir_deets) 

write.csv(PESUResults, file = file.path(wd$bin, "PESUResults2.csv"))
