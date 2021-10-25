
if(!exists("mydata_xyz")) {
  mydata_xyz <- readRDS( file.path(wd$bin, "mydata_xyz.rds") ) 
}

# Load WNS data
wns_sf <- wd$data %>%
  list.files(recursive = T, full.names = T, pattern = "shp") %>% 
  grep(pattern = "WNS_Status_Provisional", value = T) %>% 
  st_read() %>% 
  st_transform(myCRS) %>% 
  filter(
    WNS_STATUS == "Confirmed",
    !( YR_CONFIRM %in% c("2017-18", "2018-19", "2019-20") )
  )

# Newest plan -------------------------------------------------------------
set.seed(42)
subsampledOrigins <- mydata_xyz %>% 
  # METHOD IS CURRENTLY NORMALIZED PROBABILITIES
  filter(Location == "FL", analysisLab == "CASIF", method == "raw") %>% 
  group_by(ID) %>% 
  sample_n(10000, replace = T, weight = value)

distFromLineToNearestPdPositiveCounty <- pbmcapply::pbmclapply(
  1:nrow( subsampledOrigins ), 
  mc.cores = parallel::detectCores(),
  function(i) {
    p1 <- cbind(subsampledOrigins$x[i], subsampledOrigins$y[i])
    p2 <- cbind(subsampledOrigins$metersLongitude[i], subsampledOrigins$metersLatitude[i])
    myline <- rbind(p1, p2) %>% 
      st_linestring() %>% 
      st_sfc(crs = myCRS)
    whichIsClosest <- st_nearest_feature(x = myline, y = wns_sf)
    theClosest <- slice(wns_sf, whichIsClosest)
    minDistToLine <- st_distance(x = myline, y = theClosest, by_element = TRUE)
    units(minDistToLine) <- units::as_units("km") # Convert to km.
    return(minDistToLine)
    }
  )
 
subsampledOrigins_withMinDistance <- data.frame(
  subsampledOrigins, 
  distFromLineToNearestPdPositiveCounty = unlist(distFromLineToNearestPdPositiveCounty)
  )

saveRDS(subsampledOrigins_withMinDistance, file = file.path(
  wd$bin, "subsampledOrigins_withMinDistance.rds"
))

subsampledOrigins_withMinDistance %>% write.csv(file = file.path(wd$figs, "pD_distanceresults.csv"))
