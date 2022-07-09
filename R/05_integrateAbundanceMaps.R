library(lubridate)
library(terra)
suppressWarnings(suppressMessages({
  source("~/RosyFinchesdD/R/01_loadIsotopeData.R")
}))
myres <- "lr"

for(i in 1:nrow(mydata)) {

  print(paste("Working on", mydata[[i,"ID"]]))

  # Make abundance surfaces -------------------------------------------------

  if( mydata[[i,"Species"]] == "BLRF") {
    ebirdCode <- "bkrfin"
  } else if( mydata[[i,"Species"]] == "GCRF") {
    ebirdCode <- "gcrfin"
  }

  targetPath <- grep(list.dirs(wd$abundance, recursive = F), pattern = ebirdCode, value = T)

  if( mydata[[i, "sampleType"]] == "Feather" ) {
    # If sample is feather, use breeding abundance.
    myAbund0 <- file.path(targetPath, "abundance_seasonal") %>%
      list.files(pattern = myres, full.names = T) %>%
      grep(pattern = "_breeding.tif", value = T) %>%
      terra::rast()

  } else if( mydata[[i, "sampleType"]] == "Claw" ) {
    # If sample is claw, find the weeks 2-5 months prior to sampling.
    window0 <- mydata[[i, "date"]] - months(c(5,2))
    if(any(is.na(window0))) {
      # Sometimes leap years / differential month lengths cause issues.
      # Tweak by finding the first day of the week.
      window0 <- floor_date(mydata[[i, "date"]], unit = "week") - months(c(5,2))
    }
    stopifnot(!any(is.na(window0)))

    startLayer <- week(window0[[1]])
    endLayer   <- week(window0[[2]])

    if(endLayer < startLayer) {
      # If the window wraps past week 52, adjust.
      targetLayers <- c(startLayer:52, 1:endLayer)
    } else {
      targetLayers <- c(startLayer:endLayer)
    }

    ms <- file.path(targetPath, "weekly_cubes") %>%
      list.files(pattern = myres, full.names = T) %>%
      grep("abundance_median", ., value = T) %>%
      terra::rast()
    stopifnot(nlyr(ms) == 52) # 52 weeks / year expected.

    myLayers <- ms[[targetLayers]]
    # # Find 95th quantile
    # myAbund0 <- terra::quantile(myLayers, probs = c(0.95))
    # Find mean
    myAbund0 <- terra::app(myLayers, mean)

  }

  # Combine with continuous assignment maps --------------------------------

  map_dD <- file.path( wd$out, "continuous_maps", paste0(mydata[[i,"ID"]], ".tif")) %>%
    terra::rast()

  # Normalize abundance surface to sum to 1.
  myAbund1 <- myAbund0 %>%
    raster::raster() %>%
    .findProductThenNormalize() %>%
    terra::rast()

  # Resample.
  myAbund <- terra::resample(myAbund1, map_dD)

  # Find product of two surfaces then odds ratio.
  combo <- {myAbund * map_dD} %>%
    raster::raster() %>%
    isocat::.findProductThenNormalize()

  combo1 <- combo %>%
    isocat::makeOddsSurfaces()

  # Binary surface combination at >0 abundance and <2/3 odds ratio:
  combo2 <- classify(rast(combo1), c(-Inf, (2/3), 0, (2/3), Inf, 1), include.lowest=TRUE, brackets=TRUE) %>%
    raster::raster()

  # Export ------------------------------------------------------------------

  outpath <- file.path(wd$out, "combined_maps")
  if(!dir.exists(outpath)) {dir.create(outpath)}

  terra::writeRaster(combo,  overwrite = T, filename = file.path(outpath, paste0(mydata[[i, "ID"]], "_cont.tif")))
  terra::writeRaster(combo1, overwrite = T, filename = file.path(outpath, paste0(mydata[[i, "ID"]], "_OR.tif")))
  terra::writeRaster(combo2, overwrite = T, filename = file.path(outpath, paste0(mydata[[i, "ID"]], "_OR2_1.tif")))


  # Fortify and save --------------------------------------------------------
  if(!dir.exists(file.path(wd$bin, "tmp2"))) dir.create(file.path(wd$bin, "tmp2"))
  mdf <- SDMetrics::surface2df(combo1)
  saveRDS(mdf, file = file.path(wd$bin, "tmp2", paste0("abund_orig_OR_",  mydata[[i,"ID"]], ".rds")))

  # Plot --------------------------------------------------------------------
  states <- readRDS(file.path(wd$bin, "states.rds"))
  mydata_sf <- mydata[i, ] %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ) %>%
    st_transform(myCRS)
  samplesSitesx <- unique(unlist( lapply(1:length(mydata_sf$geometry), function(i) mydata_sf$geometry[[i]][1]) ))
  sampleSitesy  <- unique( unlist( lapply(1:length(mydata_sf$geometry), function(i) mydata_sf$geometry[[i]][2]) ) )

  # mdf2 <- mdf %>%
  #   dplyr::mutate(value = factor(case_when(value < 1/3 ~ 0, TRUE ~ 1)))

  p <- ggplot() +
    geom_sf( states, mapping = aes(), fill = "grey99", size = 0.25) +
    geom_tile(mdf, mapping = aes(color = value, fill = value, x=x, y=y)) +
    scale_fill_viridis_c( option = "turbo") +
    scale_color_viridis_c(option = "turbo") +
    # scale_fill_manual( breaks = c(0,1), labels = c(0,1), values = c("grey50", "goldenrod")) +
    # scale_color_manual(breaks = c(0,1), labels = c(0,1), values = c("grey50", "goldenrod")) +
    geom_sf(states, mapping = aes(), fill = NA, size = 0.25) +
    geom_star(
      data = data.frame( x=samplesSitesx, y=sampleSitesy ),
      mapping = aes(x=x,y=y),
      fill = "white", size = 3
    ) +
    theme_minimal() +
    xlab(NULL) +
    ylab(NULL) +
    theme(
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.title = element_text(size = 10)
    ) +
    theme(
      #panel.background = element_rect(fill = "grey90"),
      axis.text = element_text(size=10),
      title = element_text(size = 13),
      legend.direction = "horizontal",
      legend.title.align = c(0.5),
      legend.key.width = unit(0.8, "cm")
    )
  ggsave(p, filename = file.path(wd$figs, paste0(mydata[[i, "ID"]], "_abund_binary.png" )), width = 6)
}


# Elevation data ----------------------------------------------------------

elev0 <- raster::raster(file.path(wd$bioclim, "wc2.1_2.5m_elev.tif")) %>%
  raster:: projectRaster(crs = myCRS) %>%
  raster::crop(my_extent_aea)
elev <- elev0 %>%
  SDMetrics::surface2df()

for(i in 1:nrow(mydata)) {
  print(paste("Working on", mydata[[i,"ID"]], " - number ", i))

  combo2 <- raster::raster( file.path(outpath, paste0(mydata[[i, "ID"]], "_OR2_1.tif")))
  poly <- combo2 %>% rasterToPolygons( dissolve=TRUE) %>% st_as_sf()
  names(poly)[1] <- "layer"

  # Find everything NOT considered as a potential origin.
  antiPoly <- st_difference(st_union( states), st_union(poly))

  myxlim <- c(extent(combo2)[1] - 2e5, extent(combo2)[2] + 3e5)
  myylim <- c(extent(combo2)[3] - 1e5, extent(combo2)[4] + 1e5)


  p <- ggplot() +
    geom_sf(
      dplyr::filter(poly, layer == 1),
      mapping = aes(fill = factor(layer)),
      fill = "goldenrod", color = "goldenrod",
      ) +
    # Darken area outside boundary of analysis.
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    # Hill shading
    geom_raster(elev, mapping = aes(alpha = value, x=x, y=y), fill = "grey20") +
    scale_alpha(range =  c(0.0, 0.80), guide = "none") +
    # Show state boundaries.
    geom_sf(states, mapping = aes(), fill = NA, size = 0.25) +
    # Highlight sample site
    geom_star(
      data = data.frame( x=samplesSitesx, y=sampleSitesy ),
      mapping = aes(x=x,y=y),
      fill = "white", size = 3, alpha = 0.5, color = "black"
    ) +
    # Extra details
    theme_minimal() +
    xlab(NULL) +
    ylab(NULL) +
    coord_sf(xlim = myxlim, ylim = myylim) +
    theme(
      legend.position = "none",
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.title = element_text(size = 10)
    ) +
    theme(
      #panel.background = element_rect(fill = "grey90"),
      axis.text = element_text(size=10),
      title = element_text(size = 13),
      legend.direction = "horizontal",
      legend.title.align = c(0.5),
      legend.key.width = unit(0.8, "cm")
    )

  ggsave(p, filename = file.path(wd$figs, paste0(mydata[[i, "ID"]], "_binary_hillshade.png" )), width = 6)

}
