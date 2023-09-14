
# Setup -------------------------------------------------------------------
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopeData.R"))
keratin_augsep    <- raster::raster(file.path(wd$bin, "keratin_augsep_isoscape.tif"))
keratin_augsep_se <- raster::raster( file.path(wd$bin, "iso_se_augsep.tif") )

# Load rangemaps for each species.
BLRF_range <- raster::raster( file.path(wd$bin, "BLRF_rangeRaster.tif") )
GCRF_range <- raster::raster( file.path(wd$bin, "GCRF_rangeRaster.tif") )

# Load model residuals.
transferFunctionResids <- readRDS(file.path(wd$bin, "transferFunctionResids.rds"))

# Specify map path.
mapPath <- file.path( wd$out, "myMaps")
if(!dir.exists(mapPath)) dir.create(mapPath)

# Specify locations of fortified rasts
wd$fortified_rasts <- file.path(wd$bin,"fortified_rasts")
if(!dir.exists(wd$fortified_rasts) ) dir.create(wd$fortified_rasts)

# Create maps. -----------------------------------------------------------

lapply(c("Feather", "Claw"), function(mySampleType){
  lapply(c("BLRF", "GCRF"), function(mySpecies){

    message(paste0(mySampleType, "_", mySpecies))
    df <- dplyr::filter(mydata, Species == mySpecies, sampleType == mySampleType)
    message(nrow(df))

    if(mySpecies == "BLRF") { myRangeMap <- BLRF_range }
    if(mySpecies == "GCRF") { myRangeMap <- GCRF_range }

    if(nrow(df) > 0) {
      #for(i in 1:nrow(df)) {
      pbmcapply::pbmclapply(1:nrow(df), mc.cores = 3, function(i) {

        myMap <- isocat::isotopeAssignmentModel(
          ID               = df[i, "ID"],
          isotopeValue     = df[i, "d2H"],
          SD_indv          = transferFunctionResids,
          precip_raster    = keratin_augsep,
          precip_SD_raster = keratin_augsep_se,
          additionalModels = myRangeMap
        )
        writeRaster(myMap, filename = file.path(mapPath, paste0(df[i, "ID"], ".tif")), overwrite = T)

        # Fortify to data frame.
        mdf <- myMap %>%
          SDMetrics::surface2df() %>%
          dplyr::mutate( ID = unlist(df[i, "ID"]) )

        fwrite(mdf, file = file.path(wd$fortified_rasts, paste0(df[i, "ID"], ".csv")), row.names = F)

      })
    }
  })
})

