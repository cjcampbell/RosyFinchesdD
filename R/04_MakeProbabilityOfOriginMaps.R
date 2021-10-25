
# Setup -------------------------------------------------------------------
source('~/PESU_migration_project_directory/PESU_migration/R/00_Setup.R')
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
load( file.path(wd$bin, "range_raster.Rdata" ) )
load( file.path(wd$bin, "bestFitIso.rdata") )

# Set to 'FALSE' is you don't want to run in parallel.
nclusters <- parallel::detectCores() 

mypath <- file.path( wd$bin, "maps")
if(!dir.exists(mypath)) dir.create(mypath)

# Create maps. -----------------------------------------------------------

isocat::isotopeAssignmentModel(
  ID               = mydata_transformed$ID,
  isotopeValue     = mydata_transformed$fur_adjusted,
  SD_indv          = mydata_transformed$sdResid,
  precip_raster    = bestFitIso$isoscape,
  precip_SD_raster = bestFitIso$sd,
  savePath         = mypath,
  additionalModels = range_raster,
  nClusters = nclusters
)

# Convert maps to data.frame format for later use -------------------------

# Load maps.
maps_cropped  <- list.files(
  mypath, pattern = "Combined.*.grd$", full.names = TRUE) %>%
  raster::stack()

# Also calculate probability quantiles.
maps_quantile_stack <- lapply(1:nlayers(maps_cropped), function(i){
  isocat::makeQuantileSurfaces(maps_cropped[[i]])
}) %>% stack()
names(maps_quantile_stack) <- paste0(names(maps_cropped), "_quantile")
writeRaster(maps_quantile_stack, filename = file.path(mypath, "quantileProbabilityMaps.grd"), overwrite = TRUE)

# And odds ratios.
maps_odds_stack <- lapply(1:nlayers(maps_cropped), function(i){
  isocat::makeOddsSurfaces(maps_cropped[[i]])
}) %>% stack()
names(maps_odds_stack) <- paste0(names(maps_cropped), "_OR")
writeRaster(maps_odds_stack, filename = file.path(mypath, "ORProbabilityMaps.grd"), overwrite = TRUE)

# Make shp files for some individual bats: ----------
if(!exists("maps_odds_stack") ) maps_odds_stack <- stack(file.path(mypath, "ORProbabilityMaps.grd"))
thoseBats <- c(
  "X121.Falli",
  "X94.Mille",
  "X166.ALA28",
  "X154.Squat",
  "X91.Mille"
)
lapply(thoseBats, function(i) {
  thisLayer <- maps_odds_stack[[ which( names( maps_odds_stack )  == paste0(i, "_OR") ) ]]
  writeRaster(thisLayer, filename = file.path(mypath, paste0(i, ".tif" )), overwrite = T)
})


# Combine, make dataframe. -----------
# Do this in batches b/c it's pretty resource-intensive.
names(maps_cropped) <- paste0(names(maps_cropped), "_raw")

wd$tmp_df <- file.path(wd$bin,"tmp_df")
if(!dir.exists(wd$tmp_df) ) dir.create(wd$tmp_df)

maps_cropped_df_list <- pbmcapply::pbmclapply(
  1:nlayers(maps_cropped), mc.cores = 4, function(i){
    mdf <- stack(maps_cropped[[i]], maps_quantile_stack[[i]], maps_odds_stack[[i]]) %>%
      # So raster::as.data.frame has given me TWO big troubles here.
      # Some weird bug in raster::data.frame is messing up column names when long = TRUE.
      # AND if na.rm = TRUE, it seems to throw out *any* cell with an NA, not just a particular cell with an NA.
      # Very unhelpful if you have different ranges in your stack!
      raster::as.data.frame(xy = TRUE, long = FALSE, na.rm = FALSE) %>%
      # Hackey fixes:
      tidyr::pivot_longer(-c("x", "y"), names_to = "layer", values_to = "value") %>%
      dplyr::filter(!is.na(value)) %>%
      tidyr::separate(col = layer, into = c("ID", "method"), sep = "_")
    saveRDS(mdf, file = file.path(wd$tmp_df, paste0("df_list_", i, ".rds")))
})

maps_df <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
  lapply(readRDS) %>%
  plyr::ldply()

# Save.
saveRDS(maps_df, file = file.path(wd$bin, "maps_df.rds"))
