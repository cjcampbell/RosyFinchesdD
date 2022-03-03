
# Setup -------------------------------------------------------------------
download_GADM <- FALSE

# This script assumes that there are candidate isoscapes downloaded from isoMAP
# somewhere in the wd$data directory.
reload_isoscapes <- FALSE

# This script assumes that there are IUCN rangemaps somewhere in the wd$iucn
# directory.
reload_ebird_abundancemaps <- FALSE

# Load GADM data ----------------------------------------------------------
if(download_GADM == TRUE){

  library(rmapshaper)
  message("Loading GADM data...")

  # First, make NoAm_sf (an "sf" object) fit extent object (my_extent).
  # my_extent_aea_st <- st_bbox(st_transform(st_as_sfc(st_bbox(my_extent, crs = 4326)), myCRS))
  # saveRDS(my_extent_aea_st, file = file.path(wd$bin, "my_extent_aea_st.rds"))

  # Get GADM data to country level.
  USA <- raster::getData('GADM', path = wd$GADM, country='USA', level=0)
  CAN <- raster::getData('GADM', path = wd$GADM, country='CAN', level=0)
  MEX <- raster::getData('GADM', path = wd$GADM, country='MEX', level=0)

  # Combine into one polygon, convert to sf object.
  NoAm0 <- raster::bind( USA, CAN, MEX ) %>%
    sf::st_as_sf(.) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5e3)

  # Remove water bodies.
  USA_2 <- raster::getData('GADM', path = wd$GADM, country='USA', level=2)
  CAN_2 <- raster::getData('GADM', path = wd$GADM, country='CAN', level=2)
  MEX_2 <- raster::getData('GADM', path = wd$GADM, country='MEX', level=2)

  waterbodies <- lapply(
    list(USA_2,MEX_2,CAN_2),
    function(x){ x[(x$ENGTYPE_2) == "Water body",] }) %>%
    do.call(rbind, .) %>%
    sf::st_as_sf() %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5000)

  NoAm <- rmapshaper::ms_erase(NoAm0, waterbodies)  # Remove water bodies
  saveRDS(NoAm, file = file.path(wd$bin, "NoAm.rds"))

  # Save bufferred NoAm boundary to constrain isoscapes, etc.
  NoAm_boundary_aea <- NoAm %>%
    st_buffer(dist = 5e4) %>%
    rmapshaper::ms_erase(., waterbodies)

  saveRDS(NoAm_boundary_aea, file = file.path(wd$bin, "NoAm_boundary_aea.rds"))

  # Also create maps of USA and CAN states for detailed plotting.
  USA_1 <- raster::getData('GADM', path = wd$GADM, country='USA', level=1)
  CAN_1 <- raster::getData('GADM', path = wd$GADM, country='CAN', level=1)
  states <- raster::bind( USA_1, CAN_1) %>%
    sf::st_as_sf(.) %>%
    st_transform(crs = myCRS) %>%
    st_crop(., my_extent_aea) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 1e3)
  saveRDS(states, file = file.path(wd$bin, "states.rds"))

} else message("Not redownloading GADM Data...")


# Load isoscapes ----------------------------------------------------------

if(reload_isoscapes == TRUE){

  message("reloading isoscapes...")

  # This code requires an utd version of isocat for the `makeMultiMonthIsoscape` function.
  # Sys.unsetenv("GITHUB_PAT")
  # devtools::install_github("cjcampbell/isocat", upgrade = "never")

  projectExtendCropRaster <- function(rast) {
   rast %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea )
  }

  makeSubsettedSurfaces <- function(monthVector, rastName) {
    rast_stack <- raster::stack()
    rast_se_stack <- raster::stack()
    rast_precip <- raster::stack()
    for(i in monthVector) {
      if(i < 10) i <- paste0("0", i)
      rast_stack <- raster::stack(rast_stack, file.path(wd$isoscapes,"GlobalPrecip", paste0("d2h_", i, ".tif")))
      rast_se_stack <- raster::stack(rast_se_stack, file.path(wd$isoscapes,"GlobalPrecip", paste0("d2h_se_", i, ".tif")))
      rast_precip <- raster::stack(rast_precip, file.path(wd$precip, paste0("wc2.1_2.5m_prec_", i, ".tif")))
      }
    rast_wgs84 <- isocat::makeMultiMonthIsoscape(rast_stack, rast_se_stack, rast_precip)
    outrast <- lapply(rast_wgs84, projectExtendCropRaster)
    writeRaster(outrast[[1]], filename = file.path(wd$bin, paste0("iso_", rastName, ".tif")), overwrite = T)
    writeRaster(outrast[[2]], filename = file.path(wd$bin, paste0("iso_se_", rastName, ".tif")), overwrite = T)
  }

  # August-September surfaces
  makeSubsettedSurfaces(monthVector = 8:9, rastName = "augsep")

  # Sep-Dec surfaces
  makeSubsettedSurfaces(monthVector = 9:12, rastName = "sepdec")

  # Oct-Jan
  makeSubsettedSurfaces(monthVector = c(10:12,1), rastName = "octjan")

  # Nov-Feb
  makeSubsettedSurfaces(monthVector = c(11:12,1:2), rastName = "novfeb")

  # ## Load mean growing season isoscape and error surfaces. Combine and save. ----
  # iso_GS <- raster::raster(file.path(wd$isoscapes, "GlobalPrecipGS", "d2h_GS.tif")) %>%
  #   projectExtendCropRaster()
  # writeRaster(iso_GS, filename = file.path(wd$bin, "iso_GS.tif"), overwrite = T)
  # iso_GS_se <- raster::raster(file.path(wd$isoscapes, "GlobalPrecipGS", "d2h_se_GS.tif"))%>%
  #   projectExtendCropRaster()
  # writeRaster(iso_GS_se, filename = file.path(wd$bin, "iso_GS_se.tif"), overwrite = T)

} else message("Not reloading isoscapes, loading saved version...")

# Load ebird abundance maps -----------------------------------------------------
if(reload_ebird_abundancemaps == TRUE){
  load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
  NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )

  # Function to load maps, convert to sf, reproject, convert to raster objects.
  rangeMapConversion <- function(speciesCode) {

    # Load eBird seasonal abundance estimates.
    if(speciesCode == "BLRF") ebirdCode <- "bkrfin"
    if(speciesCode == "GCRF") ebirdCode <- "gcrfin"
    abundancePath <- file.path(wd$data, paste0(ebirdCode, "_breedingSeasonalAbundance.tif"))
    if(!file.exists(abundancePath)) {
      speciesCode_path <- ebirdst_download(species = ebirdCode, tifs_only = TRUE, force= FALSE)
      speciesCode_abund <- load_raster("abundance_seasonal", path = speciesCode_path)
      speciesCode_abund2 <- raster::crop(speciesCode_abund, extent(c(-20015109, 597136.6, 0e7, 10007555)))       # Crop to approximate area of northwestern North America for quicker reprojection
      speciesCode_aea1 <- projectRaster(speciesCode_abund2[[1]], crs = myCRS)
      speciesCode_aea <- raster::crop(speciesCode_aea1, my_extent_aea)
      writeRaster(speciesCode_aea, filename = file.path(wd$data, paste0(ebirdCode, "_breedingSeasonalAbundance.tif")), overwrite = T)
    }

    # Convert areas with non-zero expected abundance to points.
    # Draw a convex hull around said points.
    # Buffer by 100km.
    r <- raster::raster(abundancePath)
    r[r==0] <- NA
    abun_df <- raster::rasterToPoints(r) %>%
      as.data.frame()
    abun_sf <- st_as_sf(abun_df, coords = c("x", "y"))
    st_crs(abun_sf) <- myCRS
    myhull <- abun_sf %>%
      st_union() %>%
      st_convex_hull() %>%
      st_buffer(200e3)
    saveRDS(myhull, file = file.path(wd$bin, paste0(speciesCode, "_bufferedConvexHull.rds")))

    # Load candidate rangemaps.
    # Convert to simple features object and reproject to myCRS.
    myRange <- list.files(wd$bin, full.names = T, recursive = T) %>%
      grep(pattern = paste0(speciesCode, "_bufferedConvexHull.rds"), value = T) %>%
      readRDS()

    # Convert buffered rangemaps to rasters with appropriate
    ex_rast <- raster::raster( file.path(wd$bin, "augsep_iso.tif"))
    ex_rast[] <- 1
    range_raster <- raster::mask(
      ex_rast,
      mask = as_Spatial(myRange),
      updatevalue = NA
    ) %>%
      raster::mask(., NoAm_boundary_aea) %>%
      raster::crop(., my_extent_aea)

    writeRaster(range_raster, filename =  file.path(wd$bin, paste0(speciesCode, "_rangeRaster.tif")), overwrite = T)

  }

  lapply(c("BLRF", "GCRF"), rangeMapConversion)

} else message("Not reloading IUCN rangemaps...")
