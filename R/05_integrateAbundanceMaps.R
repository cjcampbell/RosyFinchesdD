library(lubridate)
library(terra)
library(sf)

suppressWarnings(suppressMessages({
  source("~/RosyFinchesdD/R/01_loadIsotopeData.R")
}))
myres <- "hr"

wd$fortified_rasts <- file.path(wd$bin,"fortified_rasts")

# Prepare seasonal abundance maps.
ex_d2H <- file.path( wd$out, "continuous_maps", paste0(myID, ".tif"))[1] %>% terra::rast()
list.dirs(wd$abundance, recursive = T) %>%
  list.files(pattern = myres, full.names = T) %>%
  grep(pattern = "bkrfin|gcrfin", value = T) %>%
  grep(pattern = "abundance_seasonal_mean", value = T) %>%
  grep(pattern = myres, value = T) %>%
  lapply(function(x) {
    terra::rast(x) %>%
      terra::project(myCRS) %>%
      terra::resample(ex_d2H, filename = file.path(wd$bin, paste0(basename(x))) )
  })

# Prepare all weekly (claw) abundance maps.
# mydata %>%
#   data.frame %>%
#   dplyr::filter(sampleType == "Claw")
# window0 <-  mydata[[ "sampleType"]] == "Claw" %>%
#   dplyr::filter() %>%
#   {.[[, "date"]] - months(c(5,2))}
# if(any(is.na(window0))) {
#   # Sometimes leap years / differential month lengths cause issues.
#   # Tweak by finding the first day of the week.
#   window0 <- floor_date(mydata[[i, "date"]], unit = "week") - months(c(5,2))
# }
# stopifnot(!any(is.na(window0)))
#
# startLayer <- week(window0[[1]])
# endLayer   <- week(window0[[2]])






# Combine abundance and d2H surfaces --------------------------------------


for(i in 1:nrow(mydata)) {

  myID <-  mydata[[i,"ID"]]

  print(paste("Working on", myID , " - number", i))

  # Make abundance surfaces -------------------------------------------------

  if( mydata[[i,"Species"]] == "BLRF") {
    ebirdCode <- "bkrfin"
  } else if( mydata[[i,"Species"]] == "GCRF") {
    ebirdCode <- "gcrfin"
  }

  if( mydata[[i, "sampleType"]] == "Feather" ) {
    next
    # If sample is feather, use breeding abundance.
    myAbund0 <- list.files(wd$bin, pattern = myres, full.names = T) %>%
      grep(pattern = ebirdCode, value = T) %>%
      grep(pattern = "abundance_seasonal_mean", value = T) %>%
      grep(pattern = myres, value = T) %>%
      terra::rast() %>%
      {.$breeding}

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

    myname <- file.path(wd$bin, paste0(ebirdCode, "_abundance_weekly_", startLayer, "_", endLayer, "_hr_2021.tif"))

    if(!file.exists(myname)) {

      if(endLayer < startLayer) {
        # If the window wraps past week 52, adjust.
        targetLayers <- c(startLayer:52, 1:endLayer)
      } else {
        targetLayers <- c(startLayer:endLayer)
      }

      ms <- file.path(targetPath, "weekly") %>%
        list.files(pattern = myres, full.names = T, recursive = T) %>%
        grep("abundance_median", ., value = T) %>%
        terra::rast()
      stopifnot(nlyr(ms) == 52) # 52 weeks / year expected.
      myLayers <- ms[[targetLayers]]
      # Find mean
      myAbund0 <- terra::app(myLayers, mean)
      map_dD <- file.path( wd$out, "continuous_maps", paste0(myID, ".tif")) %>%
        terra::rast()
      # Resample.
      myAbund <- myAbund0 %>%
        terra::project(myCRS) %>%
        terra::resample(map_dD, filename = myname)
    } else {
      myAbund <- terra::rast(myname)
    }
  }

  # Transform.
  ## Remove 0's.
  myAbund[myAbund == 0 ] <- NA

  # Fortify.
  mdf <- c(map_dD, myAbund)  %>%
    as.data.frame(xy=T)
  names(mdf)[3:4] <- c("value", "abund")

  # Fortify then calculate statistics...
  logisticCurveIt <- function(x, L, k, x0) { L / (1+exp(-k*(x-x0))) }
  combo1 <- mdf %>%
    dplyr::mutate(
      abund = abund/sum(abund, na.rm = T), # Normalize surface.
      abund_smoothed = logisticCurveIt(abund,  L = 1, k = 8, x0 = 0.25), # Transform using logistic curve (this lowers values below the 25th percentile and increases those above to close to 1.)
      combo1 = value * abund_smoothed, # Find product.
      combo1 = combo1/sum(combo1, na.rm = T) # Normalize.
    ) %>%
    arrange(abund_smoothed) %>%
    dplyr::mutate(
      combo2_csum = cumsum(combo1),
      combo2_csum_80 = case_when(combo2_csum >= 0.80 ~ 1, TRUE ~ 0),
      combo2_OR = (combo1/(1 - combo1))/(max(combo1, na.rm = T)/(1 - max(combo1, na.rm = T))),
      combo2_OR21 = case_when(combo2_OR >= 2/3 ~ 1, TRUE ~ 0),
    )

  # Export ------------------------------------------------------------------

  outpath <- file.path(wd$bin, "combined_maps")
  if(!dir.exists(outpath)) {dir.create(outpath)}

  fwrite(combo1, file = file.path(outpath, paste0(myID, ".csv")), row.names = F)

}

# New plots ----------------------------------------------------------------

states <- readRDS(file.path(wd$bin, "states.rds"))
elev0 <- list.files(wd$bioclim, pattern = "elev.tif", full.names = T, recursive = T) %>%
  terra::rast() %>%
  terra::project(myCRS) %>%
  terra::crop(my_extent_aea) %>%
  raster::raster()
elev <- elev0 %>%
  SDMetrics::surface2df()

states <- readRDS(file.path(wd$bin, "states.rds"))
hull_BLRF <- readRDS(file.path(wd$bin, paste0("BLRF", "_bufferedConvexHull.rds")))
antiPoly_BLRF <- st_difference(st_union( states), hull_BLRF)
hull_GCRF <- readRDS(file.path(wd$bin, paste0("GCRF", "_bufferedConvexHull.rds")))
antiPoly_GCRF <- st_difference(st_union( states), hull_GCRF)

plotDeets <- list(
  # Hill shading
  geom_raster(elev, mapping = aes(alpha = value, x=x, y=y), fill = "grey20"),
  scale_alpha(range =  c(0.0, 0.80), guide = "none"),
  # Show state boundaries.
  geom_sf(states, mapping = aes(), fill = NA, size = 0.25),
  guides(
    color = guide_colorbar(title.position = "top"),
    fill  = guide_colorbar(title.position = "top")
  ) ,
  theme_minimal(),
  xlab(NULL),
  ylab(NULL),
  theme(
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    axis.text = element_text(size=10),
    title = element_text(size = 13),
    legend.direction = "horizontal",
    legend.title.align = c(0.5)
  )
)

for(i in 1:nrow(mydata)) {

  myID <-  mydata[[i,"ID"]]
  print(paste("Working on",myID,"-", i))

  mdf <- file.path(wd$bin, "combined_maps", paste0(myID, ".csv")) %>%
    fread() %>%
    dplyr::filter(!is.na(value))
  mydata_sf <- mydata[i, ] %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ) %>%
    st_transform(myCRS)

  if(mydata[[i,"Species"]] == "BLRF") {
    moreDeets <- list(
      coord_sf( xlim = c(-0.5e6, 1.4e6), ylim = c(-1.4e6, 0.3e6)) ,
      theme(legend.position = c(0.85,0.90))
    )
    antiPoly <- antiPoly_BLRF
  } else {
    moreDeets <- list(
      coord_sf( xlim = c(-2.6e6, 1.5e6), ylim = c(-1.4e6, 3e6)),
      theme(legend.position = c(0.15,0.12))
    )
    antiPoly <- antiPoly_GCRF
  }

  samplesSites <- cbind(mydata[[i,"decimalLongitude"]]  , mydata[[i,"decimalLatitude"]]) %>%
    sf_project(from ="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", to = myCRS)

  # Plan probability (no combo)
  p1 <- ggplot() +
    geom_tile(mdf, mapping = aes(color = value, fill = value, x=x, y=y)) +
    scale_fill_viridis_c(  "Probability of origin", option = "turbo", na.value = NA) +
    scale_color_viridis_c( "Probability of origin", option = "turbo", na.value = NA)  +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    # Highlight sample site
    geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
    moreDeets
  ggsave(p1, file = file.path(wd$figs,  paste0(myID, "_probability_noAbund.png" ) ))

  # Combo - probability
  p2 <- ggplot() +
    geom_tile(mdf, mapping = aes(color = combo1, fill = combo1, x=x, y=y)) +
    scale_fill_viridis_c(  "Probability of origin", option = "turbo", na.value = NA) +
    scale_color_viridis_c( "Probability of origin", option = "turbo", na.value = NA)  +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    # Highlight sample site
    geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
    moreDeets
  ggsave(p2, file = file.path(wd$figs,  paste0(myID, "_probability_withAbund.png" ) ))

  # # Combo - csum
  # p3 <- ggplot() +
  #   geom_tile(mdf, mapping = aes(color = combo2_csum, fill = combo2_csum, x=x, y=y)) +
  #   scale_fill_viridis_c(  "Cumulative probability of origin", option = "turbo", na.value = NA) +
  #   scale_color_viridis_c( "Cumulative probability of origin", option = "turbo", na.value = NA)  +
  #   geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
  #   plotDeets +
  #   # Highlight sample site
  #   geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
  #   moreDeets
  # ggsave(p3, file = file.path(wd$figs,  paste0(myID, "_CumulativeProbability_withAbund.png" ) ))
  #
  # # Combo - csum binary
  # p4 <- ggplot() +
  #   geom_tile(mdf, mapping = aes(color = combo2_csum_80, fill = combo2_csum_80, x=x, y=y)) +
  #   scale_fill_viridis_c(  "Cumulative probability of origin", option = "turbo", na.value = NA) +
  #   scale_color_viridis_c( "Cumulative probability of origin", option = "turbo", na.value = NA)  +
  #   geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
  #   plotDeets +
  #   # Highlight sample site
  #   geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
  #   moreDeets
  # ggsave(p4, file = file.path(wd$figs,  paste0(myID, "_CumulativeProbability_binary_withAbund.png" ) ))

  # Odds Ratio
  p5 <- ggplot() +
    geom_tile(mdf, mapping = aes(color = combo2_OR, fill = combo2_OR, x=x, y=y)) +
    scale_fill_viridis_c( limits = c(0,1), "Odds of origin", option = "turbo", na.value = NA) +
    scale_color_viridis_c(limits = c(0,1), "Odds of origin", option = "turbo", na.value = NA)  +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
    moreDeets
  ggsave(p5, file = file.path(wd$figs,  paste0(myID, "_combo_OR.png" ) ))

  # Odds Ratio - Binary
  p6 <- ggplot() +
    geom_tile(dplyr::filter(mdf, combo2_OR21 == 1), mapping = aes(color = factor(combo2_OR21), fill = factor(combo2_OR21), x=x, y=y)) +
    scale_fill_manual(  breaks = 1, values = c("#7A0403")) +
    scale_color_manual( breaks = 1, values = c("#7A0403"))  +
    # scale_fill_viridis_d( "Odds of origin", option = "turbo", na.value = NA, begin = 0.2, end = 0.85) +
    # scale_color_viridis_d("Odds of origin", option = "turbo", na.value = NA, begin = 0.2, end = 0.85)  +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    geom_star( data = data.frame( x=samplesSites[1], y=samplesSites[2] ),  mapping = aes(x=x,y=y), fill = "white", size = 3, alpha = 0.75, color = "black") +
    moreDeets
  ggsave(p6, file = file.path(wd$figs,  paste0(myID, "_combo_OR_binary.png" ) ))
}

