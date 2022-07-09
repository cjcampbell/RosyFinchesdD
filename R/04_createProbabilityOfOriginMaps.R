
# Setup -------------------------------------------------------------------
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopeData.R"))
keratin_augsep    <- raster::raster(file.path(wd$bin, "keratin_augsep_isoscape.tif"))
keratin_augsep_se <- raster::raster( file.path(wd$bin, "iso_se_augsep.tif") )
keratin_sepdec    <- raster::raster(file.path(wd$bin, "keratin_sepdec_isoscape.tif"))
keratin_sepdec_se <- raster::raster( file.path(wd$bin, "iso_se_octjan.tif") )
keratin_octjan    <- raster::raster(file.path(wd$bin, "keratin_octjan_isoscape.tif"))
keratin_octjan_se <- raster::raster( file.path(wd$bin, "iso_se_octjan.tif") )
keratin_novfeb    <- raster::raster(file.path(wd$bin, "keratin_novfeb_isoscape.tif"))
keratin_novfeb_se <- raster::raster( file.path(wd$bin, "iso_se_novfeb.tif") )

claw_se  <- raster::raster( file.path(wd$bin, "iso_GS_se.tif"))

# Load rangemaps for each species.
BLRF_range <- raster::raster( file.path(wd$bin, "BLRF_rangeRaster.tif") )
GCRF_range <- raster::raster( file.path(wd$bin, "GCRF_rangeRaster.tif") )

# Load model residuals.
transferFunctionResids <- readRDS(file.path(wd$bin, "transferFunctionResids.rds"))

# Specify map path.
mapPath <- file.path( wd$out, "continuous_maps")
if(!dir.exists(mapPath)) dir.create(mapPath)

# Create maps. -----------------------------------------------------------

lapply(c("Feather", "Claw"), function(mySampleType){
  lapply(c("BLRF", "GCRF"), function(mySpecies){

    message(paste0(mySampleType, "_", mySpecies))
    df <- dplyr::filter(mydata, Species == mySpecies, sampleType == mySampleType)
    message(nrow(df))

    if(mySpecies == "BLRF") { myRangeMap <- BLRF_range }
    if(mySpecies == "GCRF") { myRangeMap <- GCRF_range }

    if(nrow(df) > 0) {
      pbmcapply::pbmclapply(1:nrow(df), mc.cores = 3, function(i) {

        if(mySampleType == "Feather") {
          myIsoscape <- keratin_augsep
          myIsoscape_sd <- keratin_augsep_se
        } else if(df[i, "date"] < "2020-03-01") {
          myIsoscape <- keratin_sepdec
          myIsoscape_sd <- keratin_sepdec_se
        } else if(df[i, "date"] < "2020-04-01") {
          myIsoscape <- keratin_octjan
          myIsoscape_sd <- keratin_octjan_se
        } else if(df[i, "date"] >= "2020-04-01") {
          myIsoscape <- keratin_novfeb
          myIsoscape_sd <- keratin_novfeb_se
        } else stop("date is not parsing correctly")

        myMap <- isocat::isotopeAssignmentModel(
          ID               = df[i, "ID"],
          isotopeValue     = df[i, "d2H"],
          SD_indv          = transferFunctionResids,
          precip_raster    = myIsoscape,
          precip_SD_raster = myIsoscape_sd,
          additionalModels = myRangeMap
        )
        writeRaster(myMap, filename = file.path(mapPath, paste0(df[i, "ID"], ".tif")), overwrite = T)
      })
    }

  })
})

# Convert maps to data.frame format for later use -------------------------

# Load maps.
maps_cropped  <- list.files(mapPath, pattern = ".tif$", full.names = TRUE) %>%
  raster::stack()
#writeRaster(maps_cropped, filename = file.path(wd$out, "normalizedProbabilityMaps.grd"), overwrite = TRUE)

writeRaster(
  maps_cropped[[unlist(mydata[mydata$Species == "BLRF", "ID"])]],
  filename = file.path(wd$out, "BLRF_normalizedProbabilityMaps.grd"),
  overwrite = T)
writeRaster(
  maps_cropped[[unlist(mydata[mydata$Species == "GCRF", "ID"])]],
  filename = file.path(wd$out, "GCRF_normalizedProbabilityMaps.grd"),
  overwrite = T)

# Make odds ratios.
maps_odds_stack <- lapply(1:nlayers(maps_cropped), function(i){
  isocat::makeOddsSurfaces(maps_cropped[[i]])
}) %>% raster::stack()
names(maps_odds_stack) <- paste0(names(maps_cropped), "_OR")
#writeRaster(maps_odds_stack, filename = file.path(wd$out, "OddsProbabilityMaps.grd"), overwrite = TRUE)
writeRaster(
  maps_odds_stack[[paste0(unlist(mydata[mydata$Species == "BLRF", "ID"]), "_OR")]],
  filename = file.path(wd$out, "BLRF_oddsProbabilityMaps.grd"),
  overwrite = T)
writeRaster(
  maps_odds_stack[[paste0(unlist(mydata[mydata$Species == "GCRF", "ID"]), "_OR")]],
  filename = file.path(wd$out, "GCRF_oddsProbabilityMaps.grd"),
  overwrite = T)


# # Combine, fortify into dataframe. -----------
# # Do this in batches b/c it's pretty resource-intensive.
# names(maps_cropped) <- paste0(names(maps_cropped), "_raw")
#
# wd$tmp_df <- file.path(wd$bin,"tmp_df")
# if(!dir.exists(wd$tmp_df) ) dir.create(wd$tmp_df)
#
# maps_cropped_df_list <- pbmcapply::pbmclapply(
#   1:nlayers(maps_cropped), mc.cores = 4, function(i){
#     mdf <- stack(maps_cropped[[i]], maps_quantile_stack[[i]], maps_odds_stack[[i]]) %>%
#       raster::as.data.frame(xy = TRUE, long = FALSE, na.rm = FALSE) %>%
#       tidyr::pivot_longer(-c("x", "y"), names_to = "layer", values_to = "value") %>%
#       dplyr::filter(!is.na(value)) %>%
#       tidyr::separate(col = layer, into = c("ID", "method"), sep = "_")
#     saveRDS(mdf, file = file.path(wd$tmp_df, paste0("df_list_", i, ".rds")))
# })
#
# maps_df <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
#   lapply(readRDS) %>%
#   plyr::ldply()
#
# # Save.
# saveRDS(maps_df, file = file.path(wd$bin, "maps_df.rds"))
#
# # Make odds-ratio 2:1 surfaces -------------------------------------------------
# ORmaps <- raster::stack(file.path(wd$bin, "ORProbabilityMaps.grd"))
# # Reclassify values below 2:1 threshold to 0, above to 1.
# ORmaps_2to1 <- raster::reclassify(ORmaps, rcl = c(-Inf,1/3,0,1/3,Inf,1))
# writeRaster(ORmaps_2to1, file = file.path(wd$out, "2-1OddsRatioSurfaces.grd"))
#
# mydata_sf <- mydata %>%
#   st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ) %>%
#   st_transform(myCRS)
# ranges <- lapply( c("BLRF", "GCRF"), function(speciesCode){
#   list.dirs(wd$iucn) %>%
#     grep(pattern = speciesCode, value = T) %>%
#     grep(pattern = "species_data", value = T) %>%
#     rgdal::readOGR(dsn = ., layer = "data_0") %>%
#     st_as_sf(crs = 4326) %>%
#     st_transform(crs = myCRS) %>%
#     st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
#     st_make_valid() %>%
#     dplyr::mutate(SEASONAL = factor(SEASONAL))
# })
#
# ## Plot each map as png-----
# if(!exists("maps_df")) maps_df <- readRDS(file.path(wd$bin, "maps_df.rds"))
# states <- readRDS(file.path(wd$bin, "states.rds"))
#
# makePlot <- function(df, myrow, continuous = TRUE) {
#
#   if(!continuous) {
#     df <- dplyr::mutate(df, value = factor(value, levels = c("0", "1"), labels = c("Low", "High")))
#   }
#
#   myxlim <- c(min(df$x) - 2e5, max(df$x) + 3e5)
#   myylim <- c(min(df$y) - 1e5, max(df$y) + 1e5)
#
#   if(myxlim[1] <= -25e5) myxlim[1] <- -25e5
#
#   p <- ggplot() +
#     geom_tile(df, mapping = aes(x=x,y=y,fill=value, color = value)) +
#     geom_star(
#       data.frame(x=myrow$geometry[[1]][1], y=myrow$geometry[[1]][2]),
#       mapping = aes(x=x,y=y),
#       fill = "white", size = 3
#     ) +
#     geom_sf(states, mapping = aes(), fill = NA, size = 0.25) +
#     theme_minimal() +
#     xlab(NULL) +
#     ylab(NULL) +
#     ggtitle(paste0(myrow$Species, " ", myrow$ID)) +
#     theme(
#       axis.text = element_text(size=10),
#       title = element_text(size = 13),
#       panel.grid = element_blank(),
#       legend.direction = "horizontal",
#       legend.title.align = c(0.5),
#       legend.key.width = unit(0.8, "cm")
#     )
#
#   if(continuous) {
#     p <- p  +
#       scale_fill_viridis_c( limits = c(0,1), "Odds of origin", option = "turbo") +
#       scale_color_viridis_c(limits = c(0,1), "Odds of origin", option = "turbo") +
#       guides(
#         color = guide_colorbar(title.position = "top"),
#         fill  = guide_colorbar(title.position = "top")
#       ) +
#       theme(
#         legend.key.height = unit(0.5, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         legend.title = element_text(size = 10)
#       )
#   } else {
#     p <- p  +
#       scale_fill_manual(  "Odds of origin" , values = c("#30123BFF", "#4BF48A"), drop = F ) +
#       scale_color_manual( "Odds of origin" , values = c("#30123BFF", "#4BF48A"), drop = F ) +
#       guides(
#         color = guide_legend(title.position = "top"),
#         fill  = guide_legend(title.position = "top")
#       ) +
#       theme(
#         text = element_text(size = 8),
#         legend.title = element_text(size = 10)
#       )
#   }
#
#   # Adjust location of legend + scale + rectangular background
#   # depending on species being plotted.
#   if(myrow$Species == "BLRF") {
#     boxmin <- myxlim[2] - ((myxlim[2]-myxlim[1]) * 0.28) + 2e4
#     boxmax <- myxlim[2] + 9e4
#     p <- p +
#       ggspatial::annotation_scale(
#         width_hint = 0.22,
#         bar_cols = c("grey50", NA),
#         location="br"
#       )
#
#     if(!continuous){
#       p <- p +  theme( legend.position = c(0.86,0.10) )
#     } else {
#       p <- p + theme(
#         legend.text = element_text(size = 8),
#         legend.position = c(0.86,0.15) )
#     }
#
#   } else if(myrow$Species == "GCRF") {
#     boxmin <- myxlim[1] - 16e4
#     boxmax <- myxlim[1] + ((myxlim[2]-myxlim[1]) * 0.28) + 15e4
#     p <- p +
#       ggspatial::annotation_scale(
#         width_hint = 0.22,
#         bar_cols = c("grey50", NA),
#         location="bl"
#       )
#     if(!continuous){
#      p <- p + theme(
#        legend.position = c(0.17,0.09)
#        )
#     } else {
#      p <- p + theme(
#        legend.position = c(0.17,0.15)
#        )
#     }
#   }
#
#   p <- p + coord_sf(xlim = myxlim, ylim = myylim)
#
#   return(p)
# }
#
# lapply(mydata$ID, function(i){
#   print(i)
#   myrow <- mydata_sf[mydata_sf$ID == i , ]
#   mydf_OR <- dplyr::filter(maps_df, ID == i, method  == "OR")
#
#   stopifnot(nrow(myrow) == 1, nrow(mydf_OR) > 1)
#
#   # Make continuous odds plots.
#   p <- makePlot(mydf_OR, myrow = myrow, continuous = T)
#   ggsave(p, file = file.path(wd$figs, paste0(myrow$ID, "_odds.png")), width = 8, height = 8)
#
#   # Make thresholded odds plots.
#   mydf_2t1 <- dplyr::mutate(mydf_OR, value = factor(case_when(value >= 1/3 ~ 1, TRUE ~ 0)))
#   p2 <- makePlot(mydf_2t1, myrow = myrow, continuous = F)
#   ggsave(p2, file = file.path(wd$figs, paste0(myrow$ID, "_2t1_odds.png")), width = 8, height= 8)
# })
#
