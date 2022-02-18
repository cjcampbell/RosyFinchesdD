
# Setup -------------------------------------------------------------------
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopeData.R"))
augsep_feather <- raster::raster(file.path(wd$bin, "featherIsoscape.tif"))
augsep_se <- raster::raster(file.path(wd$bin, "augsep_se.tif"))

# Load rangemaps for each species.
BLRF_range <- raster::raster( file.path(wd$bin, "BLRF_rangeRaster.tif") )
GCRF_range <- raster::raster( file.path(wd$bin, "GCRF_rangeRaster.tif") )

# Load model residuals.
transferFunctionResids <- readRDS(file.path(wd$bin, "transferFunctionResids.rds"))

# Specify map path.
mapPath <- file.path( wd$bin, "maps")
if(!dir.exists(mapPath)) dir.create(mapPath)

# Create maps. -----------------------------------------------------------

lapply(c("Feather", "Claw"), function(mySampleType){
  lapply(c("BLRF", "GCRF"), function(mySpecies){

    message(paste0(mySampleType, "_", mySpecies))
    df <- dplyr::filter(mydata, Species == mySpecies, sampleType == mySampleType)
    message(nrow(df))

    if(mySampleType == "Feather") {
      myIsoscape <- augsep_feather
      myIsoscape_sd <- augsep_se
    } else {
      myIsoscape <- augsep_feather
      myIsoscape_sd <- augsep_se
    }
    if(mySpecies == "BLRF") {
      myRangeMap <- BLRF_range
      myIsoscape <- myIsoscape
    } else if(mySpecies == "GCRF") {
      myRangeMap <- GCRF_range
      myIsoscape <- myIsoscape
    }

    if(nrow(df) > 0) {
      pbmcapply::pbmclapply(1:nrow(df), mc.cores = 3, function(i) {

        myMap <- isocat::isotopeAssignmentModel(
          ID               = df[i, "ID"],
          isotopeValue     = df[i, "d2H"],
          SD_indv          = transferFunctionResids,
          precip_raster    = myIsoscape,
          precip_SD_raster = myIsoscape_sd,
          additionalModels = myRangeMap,
          additionalModel_name = "eBird"
        )
        writeRaster(myMap, filename = file.path(mapPath, paste0(df[i, "ID"], ".tif")), overwrite = T)
      })
    }

  })
})


# Stack and process maps. ----

myStack <- raster::stack(list.files(mapPath, full.names = T))

mydata_sf <- mydata %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ) %>%
  st_transform(myCRS)
ranges <- lapply( c("BLRF", "GCRF"), function(speciesCode){
  list.dirs(wd$iucn) %>%
    grep(pattern = speciesCode, value = T) %>%
    grep(pattern = "species_data", value = T) %>%
    rgdal::readOGR(dsn = ., layer = "data_0") %>%
    st_as_sf(crs = 4326) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
    st_make_valid() %>%
    dplyr::mutate(SEASONAL = factor(SEASONAL))
})

if(!exists("states")) states <- readRDS(file.path(wd$bin, "states.rds"))

for(i in 1:nlayers(myStack) ) {
  print(paste0(i, " of ", nlayers(myStack)))
  myrow <- mydata_sf[mydata_sf$ID == names(myStack[[i]]) , ]
  if(!is.na(myrow$d2H)) {
    odds <- isocat::makeOddsSurfaces(myStack[[i]])

    mySpecies <- myrow$Species
    if(mySpecies == "BLRF") {
      myRangeMap <- ranges[[1]]
    } else if(mySpecies == "GCRF") {
      myRangeMap <- ranges[[2]]
    }

    myTissue <- myrow$sampleType
    if(myTissue == "Claw") {
      highlight <- NA
      lowlight <- NA
    # } else if(mySpecies == "GCRF") {
    #   highlight <- myRangeMap %>% dplyr::filter(SEASONAL %in% c(1,2))
    #   lowlight <-  myRangeMap %>%  dplyr::filter(SEASONAL %in% c(3,4))
    } else {
      highlight <- NA
      lowlight <- NA
    }

    myxlim <- c(min(df$x) - 2e5, max(df$x) + 3e5)
    myylim <- c(min(df$y) - 1e5, max(df$y) + 1e5)
    boxmin <- myxlim[2] - ((myxlim[2]-myxlim[1]) * 0.28)
    boxmax <- myxlim[2] +7e4

    df <- SDMetrics::surface2df(odds)
    p <- ggplot() +
      geom_tile(df, mapping = aes(x=x,y=y,fill=value, color = value)) +
      scale_fill_viridis_c( limits = c(0,1), "Odds of origin", option = "turbo") +
      scale_color_viridis_c(limits = c(0,1), "Odds of origin", option = "turbo") +
      geom_star(
        data.frame(x=myrow$geometry[[1]][1], y=myrow$geometry[[1]][2]),
        mapping = aes(x=x,y=y),
        fill = "white", size = 3
        ) +
      geom_sf(states, mapping = aes(), fill = NA, size = 0.25) +
      theme_minimal() +
      xlab(NULL) +
      ylab(NULL) +
      ggtitle(paste0(myrow$Species, " ", myrow$ID)) +
      theme(
        panel.grid = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.86,0.11),
        legend.title.align = c(0.5),
        legend.key.width = unit(0.8, "cm")
      ) +
      guides(
        color = guide_colorbar(title.position = "top"),
        fill  = guide_colorbar(title.position = "top")
        )  +
      geom_rect(
        fill = "white", color = "grey20",
        data = data.frame(
          xmin = boxmin,
          xmax = boxmax,
          ymin = myylim[1]-6e4,
          ymax = myylim[1] + 2e5
        ),
        mapping = aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
      ) +
      ggspatial::annotation_scale(
        width_hint = 0.22,
        bar_cols = c("grey50", NA),
        location="br"
      )

    if(!is.na(highlight)) {
      p <- p + geom_sf(highlight, mapping = aes(), fill = NA, color = "black", size = 1)
    }
    if(!is.na(lowlight)) {
      p <- p + geom_sf(lowlight, mapping = aes(), fill = "grey50", alpha = 0.5, color = NA, size = 1)
    }
    p <- p + coord_sf(xlim = myxlim, ylim = myylim)

    ggsave(p, file = file.path(wd$out, paste0(myrow$ID, ".png")), width = 8)
  }
}


# Convert maps to data.frame format for later use -------------------------

# Load maps.
maps_cropped  <- list.files(mapPath, pattern = ".tif", full.names = TRUE) %>%
  raster::stack()

# Also calculate probability quantiles.
maps_quantile_stack <- lapply(1:nlayers(maps_cropped), function(i){
  isocat::makeQuantileSurfaces(maps_cropped[[i]])
}) %>%
  raster::stack()
names(maps_quantile_stack) <- paste0(names(maps_cropped), "_quantile")
writeRaster(maps_quantile_stack, filename = file.path(wd$bin, "quantileProbabilityMaps.grd"), overwrite = TRUE)

# And odds ratios.
maps_odds_stack <- lapply(1:nlayers(maps_cropped), function(i){
  isocat::makeOddsSurfaces(maps_cropped[[i]])
}) %>% raster::stack()
names(maps_odds_stack) <- paste0(names(maps_cropped), "_OR")
writeRaster(maps_odds_stack, filename = file.path(wd$bin, "ORProbabilityMaps.grd"), overwrite = TRUE)


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
