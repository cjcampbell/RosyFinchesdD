
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

  # Function to extend/crop/mask by above.
  ECM <- function(rasterLayer){
    rasterLayer %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea ) %>%
      raster::resample(., refIsoscape) %>%
      raster::mask(   ., NoAm_boundary_aea  )
  }

  # Load then write assignR's isoscape. Could call directly, but this works okay.
  assignRpathPattern <- "assignR_d2h_world"
  # It's called something different in v. 2... grr! Will use original version already stored locally.

  # GGS_H <- assignR::d2h_lrNA
  # if( !dir.exists(file.path(wd$data, assignRpathPattern)) ) {
  #   dir.create(file.path(wd$data, assignRpathPattern))
  # }
  # writeRaster(GGS_H[[1]], overwrite = TRUE,
  #             filename = file.path(wd$data, assignRpathPattern, "isoscape.tif"))
  # writeRaster(GGS_H[[2]], overwrite = TRUE,
  #             filename = file.path(wd$data, assignRpathPattern, "sd.tif"))

  # Load then crop and mask environmental data rasters.
  loadAdjustIsoscapeTiff <- function( directory, path_pattern, isoscape_pattern,
                                      sd_pattern, refIsoscape){
    l <- list()
    l$directory          <- directory
    l$path_pattern       <- path_pattern
    l$isoscape_pattern   <- isoscape_pattern
    l$sd_pattern         <- sd_pattern
    l$reference_isoscape_name <- names(refIsoscape)

    l$isoscape <- list.files(
      directory, recursive = TRUE, pattern = isoscape_pattern, full.names = TRUE
    ) %>%
      grep(path_pattern, ., value = TRUE) %>%
      raster::raster(.) %>%
      ECM(.)

    l$sd <- list.files(
      directory, recursive = TRUE, pattern = sd_pattern, full.names = TRUE
    ) %>%
      grep(path_pattern, ., value = TRUE) %>%
      raster::raster(.) %>%
      ECM(.)

    return(l)
  }

  suppressWarnings({

    refIsoscape <- file.path(wd$data, assignRpathPattern, "isoscape.tif") %>%
      raster::raster(.) %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea )

    my_isoscapes <- mapply(
      FUN = loadAdjustIsoscapeTiff,
      path_pattern     = c("66100", "83507", assignRpathPattern),
      isoscape_pattern = c(rep("predkrig.tiff$", 2) , "isoscape.tif$"),
      sd_pattern       = c(rep("stdkrig.tiff$" , 2) , "sd.tif$"),
      MoreArgs = list(
        directory = wd$data,
        refIsoscape = refIsoscape
        ),
      SIMPLIFY = FALSE
    )

  })

  # Check that everything turned out okay.
  lapply(my_isoscapes, function(i) c( i$isoscape, i$sd) ) %>%
    unlist %>%
    lapply(., compareRaster, x = .[[1]] ) %>%
    unlist %>%
    all %>%
    {if(.!=TRUE) stop("Something is wrong! CompareRaster yeilds FALSE results.")}

  # Save.
  save(my_isoscapes, file = file.path(wd$bin, "my_isoscapes.RData"))
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
