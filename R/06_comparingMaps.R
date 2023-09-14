source("~/RosyFinchesdD/.Rprofile")

# Load abundance-weighted origin maps.
myfiles <- list.files(file.path(wd$out, "combined_maps"), recursive = TRUE, pattern = "_cont.tif$", full.names = T)
maps <- myfiles %>%
  bettermc::mclapply(.,mc.cores = 3, function(x) {
    .findProductThenNormalize( makeQuantileSurfaces(raster(x)) )
  }) %>%
  raster::stack()
names(maps) <- gsub("_cont.tif", "", basename(myfiles))

# Load data.
suppressWarnings(suppressMessages({ source("~/RosyFinchesdD/R/01_loadIsotopeData.R") }))

# Make simmatrix ------------------------------------------------
# For each species:
lapply(c("GCRF", "BLRF"), function(spp) {
  lapply(c("Claw", "Feather"), function(st) {

    whichmaps <- unlist(mydata[ mydata$Species == spp & mydata$sampleType == st, "ID" ])
    if (length(whichmaps) != 0) {
      print(paste("Working on", spp, st))
      mymaps <- maps[[ whichmaps ]]
      simmatrix <- isocat::simmatrixMaker(
        mymaps,
        nClusters = FALSE,
        csvSavePath = FALSE
      )
      saveRDS(
        simmatrix,
        file = file.path( wd$bin, paste(spp, st, "simmatrix.rds", sep ="_") )
        )
    }
  })
})

# Perform hierarchical clustering -----------------------------------------
# For each species, perform normal clustering based on similarity matrices.

lapply(c("GCRF", "BLRF"), function(spp) {
 lapply(c("Claw", "Feather"), function(st) {
   mypath <- file.path( wd$bin, paste(spp, st, "simmatrix.rds", sep = "_"))
   if(file.exists(mypath)) {
    simmatrix <- readRDS(mypath)
    simmatrix_dist <- dist( 1 - simmatrix, diag = TRUE, upper = TRUE )
    clustered_simmatrix <- hclust( simmatrix_dist, method = "ward.D2")
    save(clustered_simmatrix, file = file.path(wd$bin, paste(spp, st, "clustered_simmatrices.Rdata", sep= "_") ))

    wss <- rep(NA, 15)
    for (i in 1:15) wss[i] <- sum(kmeans(simmatrix_dist, i)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

    if(spp == "GCRF") { k <- 4} else {k <- 4}
    # Plot trees.
    plot(clustered_simmatrix, cex = 0.6)
    rect.hclust(clustered_simmatrix, k = k, border = 2:8)


    mytree <- clustered_simmatrix
    rawCluster <- stats::cutree( mytree, k = k )

    myClustDeets <- data.frame(
      rawCluster = rawCluster, ID = mytree$labels,
      Species = spp, sampleType = st
      )

    # Reorder clusters w/ respect to d2H vals (approximately, latitude)
    myClustDeets <- myClustDeets %>%
      left_join(mydata) %>%
      group_by(rawCluster) %>%
      dplyr::summarise(mean = mean(d2H)) %>%
      arrange(mean) %>%
      dplyr::mutate(cluster = row_number()) %>%
      dplyr::select(rawCluster, cluster) %>%
      left_join(myClustDeets, .) %>%
      dplyr::select(-rawCluster)

    saveRDS(myClustDeets, file = file.path(wd$bin, paste(spp, st, "clustDeets.rds", sep = "_")))
   }
 })
})

mydata_clustered <- list.files(wd$bin, pattern = "clustDeets.rds$", full.names = T) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  full_join(mydata, by = c("ID", "Species", "sampleType"))

saveRDS(mydata_clustered, file = file.path(wd$bin, "mydata_clustered.rds"))


# Exploratory analyses ----------------------------------------------------

mydata_clustered %>%
  ggplot() +
  aes(x = Species, y = cluster) +
  geom_jitter()




