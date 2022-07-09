
# Load abundance-weighted origin maps.
maps <- list.files(file.path(wd$out, "combined_maps"), recursive = TRUE, pattern = "_cont.tif$", full.names = T) %>%
  lapply(terra::rast) %>%
  terra::rast() %>%
  raster::stack()
names(maps) <- gsub("_cont", "", names(maps))

# Load data.
suppressWarnings(suppressMessages({ source("~/RosyFinchesdD/R/01_loadIsotopeData.R") }))

# Make simmatrix ------------------------------------------------
# For each species:
lapply(c("GCRF", "BLRF"), function(spp) {
  print(paste("Working on", spp))
  mymaps <- maps[[ unlist(mydata[ mydata$Species == spp, "ID" ]) ]]
  simmatrix <- isocat::simmatrixMaker(
    mymaps,
    nClusters = FALSE,
    csvSavePath = FALSE
  )
  saveRDS(simmatrix, file = file.path( wd$bin, paste0(spp, "_simmatrix.rds") ))
})

# Perform hierarchical clustering -----------------------------------------
# For each species, perform normal clustering based on similarity matrices.

lapply(c("GCRF", "BLRF"), function(spp) {
  simmatrix <- readRDS(file.path( wd$bin, paste0(spp, "_simmatrix.rds")))
  simmatrix_dist <- dist( 1 - simmatrix, diag = TRUE, upper = TRUE )
  clustered_simmatrix <- hclust( simmatrix_dist, method = "ward.D2")
  save(clustered_simmatrix, file = file.path(wd$bin, paste0(spp, "clustered_simmatrices.Rdata") ))

  wss <- rep(NA, 15)
  for (i in 1:15) wss[i] <- sum(kmeans(simmatrix_dist, i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

  if(spp == "GCRF") { k <- 3} else {k <- 2}
  # Plot trees.
  plot(clustered_simmatrix, cex = 0.6)
  rect.hclust(clustered_simmatrix, k = k, border = 2:8)

})

#  Prepare to cut trees into groups of individuals with common origins ----------------

lapply(c("GCRF", "BLRF"), function(spp) {
# k-means cluster number
wss <- rep(NA, 15)
for (i in 1:15) wss[i] <- sum(kmeans(simmatrix_dist, i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# Plot trees.
plot(clustered_simmatrix, cex = 0.6)
rect.hclust(clustered_simmatrix, k = 6, border = 2:8)


# # # Estimate the latitude of highest probability of origin for each surface ------
# mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
# maps_cropped_df    <- readRDS( file.path(wd$bin, "maps_df.rds") ) %>%
#   dplyr::filter(method == "raw")
#
# ave_lat <- lapply(1:nrow(mydata_transformed), function(i){
#   ex <- mydata_transformed[i, ]
#   u1 <- na.omit( maps_cropped_df[ maps_cropped_df$ID == ex$SampleName , ] )
#   lats <- base::sample(size = 1000, x = u1$y, replace = TRUE, prob = u1$value)
#   mean(lats)
# } ) %>%
#   unlist
#
# mydata_aveLat <- mydata_transformed %>% dplyr::mutate(ave_lat = ave_lat)
#
# saveRDS(mydata_aveLat, file = file.path(wd$bin, "mydata_aveLat.rds"))

# Manually decide how many clusters to cut into. -----------------------
# This part requires manual input!
my_k <- 4

# Cut trees into groups ---------------------------------------------------
# Load data.
mydata_aveLat <- readRDS( file.path(wd$bin, "mydata_aveLat.rds") )
load(file.path(wd$bin, "clustered_simmatrices.Rdata"))

# Apply.

mytree <- clustered_simmatrix
rawCluster <- stats::cutree( mytree, k = my_k )

myClustDeets <- as.data.frame(rawCluster)
myClustDeets$SampleName <- mytree$labels

# Reorder cluster numbers with respect to descending latitude.
myClustDeets0 <- mydata_aveLat %>%
  right_join(myClustDeets, by = "SampleName")
lat_key <- myClustDeets0 %>%
  group_by(rawCluster) %>%
  dplyr::summarise(meanLat = mean(ave_lat, na.rm = TRUE)) %>%
  arrange(meanLat)
lat_key$order <- 1:nrow(lat_key)

save(lat_key, file = file.path( wd$bin, "recodedClusters.Rdata") )

myClustDeets <- myClustDeets %>%
  dplyr::mutate(OriginCluster = plyr::mapvalues(
    rawCluster, lat_key$rawCluster, lat_key$order)
  ) %>%
  dplyr::select(SampleName, OriginCluster)

myClustDeets$originColor <- viridisLite::viridis(
  length(unique(myClustDeets$OriginCluster)))[myClustDeets$OriginCluster]

mydata_clustered <- myClustDeets

save(mydata_clustered, file = file.path(wd$bin, "mydata_clustered.Rdata"))

