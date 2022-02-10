
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
  USA_2 <- raster::getData('GADM', path = wd$bin, country='USA', level=2)
  CAN_2 <- raster::getData('GADM', path = wd$bin, country='CAN', level=2)
  MEX_2 <- raster::getData('GADM', path = wd$bin, country='MEX', level=2)

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

} else message("Not redownloading GADM Data...")


# Load isoscapes ----------------------------------------------------------

if(reload_isoscapes == TRUE){
  message("reloading isoscapes...")

  NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )

  augsep_stack <- raster::stack(
    file.path(wd$isoscapes,"GlobalPrecip", "d2h_08.tif"),
    file.path(wd$isoscapes,"GlobalPrecip", "d2h_09.tif")
  )
  augsep_ave <- raster::calc(augsep_stack, fun = mean)
  augsep_iso <- augsep_ave %>%
    raster::projectRaster(., crs = myCRS) %>%
    raster::extend( ., my_extent_aea ) %>%
    raster::crop(   ., my_extent_aea )
  writeRaster(augsep_iso, filename = file.path(wd$bin, "augsep_iso.tif"), overwrite = T)

  augsep_se_stack <- raster::stack(
    file.path(wd$isoscapes,"GlobalPrecip", "d2h_se_08.tif"),
    file.path(wd$isoscapes,"GlobalPrecip", "d2h_se_09.tif")
  )
  sqrt_mse <- function(x,y) { sqrt((x^2) + (y^2)) }
  augsep_se <- augsep_se_stack %>%
    raster::projectRaster(., crs = myCRS) %>%
    raster::extend( ., my_extent_aea ) %>%
    raster::crop(   ., my_extent_aea ) %>%
    # Calculate combined se.
    raster::overlay(., fun = sqrt_mse)
  writeRaster(augsep_se, filename = file.path(wd$bin, "augsep_se.tif"), overwrite = T)


} else message("Not reloading isoscapes, loading saved version...")

# Load ebird abundance maps -----------------------------------------------------
if(reload_ebird_abundancemaps == TRUE){
  load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
  NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )

  # Function to load maps, convert to sf, reproject, convert to raster objects.
  rangeMapConversion <- function(speciesCode) {

    # Load eBird seasonal abundance estimates.
    speciesCode_path <- ebirdst_download(species = speciesCode, tifs_only = TRUE, force= FALSE)
    speciesCode_abund <- load_raster("abundance_seasonal", path = speciesCode_path)
    speciesCode_abund2 <- raster::crop(speciesCode_abund, extent(c(-20015109, 597136.6, 0e7, 10007555)))       # Crop to approximate area of northwestern North America for quicker reprojection
    speciesCode_aea1 <- projectRaster(speciesCode_abund2[[1]], crs = myCRS)
    speciesCode_aea <- raster::crop(speciesCode_aea1, my_extent_aea)
    writeRaster(speciesCode_aea, filename = file.path(wd$data, paste0(speciesCode, "_breedingSeasonalAbundance.tif")))

    # Convert areas with non-zero expected abundance to points.
    # Draw a convex hull around said points.
    # Buffer by 100km.
    r <- raster::raster(file.path(wd$data, paste0(speciesCode, "_breedingSeasonalAbundance.tif")))
    r[r==0] <- NA
    abun_df <- raster::rasterToPoints(r) %>%
      as.data.frame()
    abun_sf <- st_as_sf(abun_df, coords = c("x", "y"))
    st_crs(abun_sf) <- myCRS
    myhull <- abun_sf %>%
      st_union() %>%
      st_convex_hull() %>%
      st_buffer(100e3)
    saveRDS(myhull, file = file.path(wd$bin, paste0(speciesCode, "_bufferedConvexHull.rds")))

    # Load candidate rangemaps.
    # Convert to simple features object and reproject to myCRS.
    myRange <- list.files(wd$bin, full.names = T, recursive = T) %>%
      grep(pattern = paste0(speciesCode, "_bufferedConvexHull.rds"), value = T) %>%
      readRDS()

    # Convert buffered rangemaps to rasters with appropriate
    ex_rast <- my_isoscapes[[1]]$isoscape
    ex_rast[] <- 1
    range_raster <- raster::mask(
      ex_rast,
      mask = as_Spatial(myRange),
      updatevalue = NA
    ) %>%
      raster::mask(., NoAm_boundary_aea) %>%
      raster::crop(., my_extent_aea)

    saveRDS(range_raster, file = file.path(wd$bin, paste0(speciesCode, "_rangeRaster.rds")))

  }

  lapply(c("BLRF", "GCRF"), rangeMapConversion)

} else message("Not reloading IUCN rangemaps...")
