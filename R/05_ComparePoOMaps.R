
# Setup --------------------------------------------------------------------

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
load( file.path(wd$bin, "range_raster.Rdata" ) )
load( file.path(wd$bin, "bestFitIso.rdata") , verbose = T)


# Make "grouped" simmatrix ------------------------------------------------
# This approach find within-species similarity values for groups of individuals
# where each group representes a distinct combination of model parameter inputs,
# i.e., individuals with identical isotope values will be considered as one
# group for the model run then the simmatrix expanded to represent them
# individually in the next step.

do_par <- TRUE
how_many_cores <- 3     # Even with this optimization,
# it's a very resource-intensive process.

# Run.
df <- mydata_transformed 

message("Creating groups for how many individuals: ", nrow(df))

df$group <- paste0( "G", group_indices(df, fur_adjusted, sdResid) )

bygroup <- df %>%
  dplyr::group_by(group) %>%
  dplyr::sample_n(1) %>%
  dplyr::ungroup()

groupModels <- isocat::isotopeAssignmentModel(
  ID               = bygroup$group,
  isotopeValue     = bygroup$fur_adjusted,
  SD_indv          = bygroup$sdResid,
  precip_raster    = bestFitIso$isoscape,
  precip_SD_raster = bestFitIso$sd,
  savePath         = FALSE,
  additionalModels = range_raster,
  nClusters = how_many_cores
)

y <- system.time(
  b <- isocat::simmatrixMaker(
    groupModels,
    nClusters = FALSE,
    csvSavePath = FALSE
  )
)
save(b, file = file.path( wd$bin, "group_simmatrix.Rdata"))
write.csv(
  cbind(individuals = nrow(df), y),
  file = file.path( wd$bin, "log.csv")
)
message("Simmatrix completed for groups." )


# Populate simmatrix ------------------------------------------------------
message("Populating simmatrix.")
load( file.path( wd$bin, "group_simmatrix.Rdata") )

# Duplicate rows when converting from groups back to individuals.
system.time({
  ids <- as.vector(df$ID)
  x <- matrix( as.numeric(NA), length(ids), length(ids), dimnames=list(ids,ids) )
  # x <- backup
  ut <- upper.tri(x, diag = TRUE)

  out <- pbmcapply::pbmclapply(1:length(x),  x = x, ut = ut, mc.cores = how_many_cores, function(i, x, ut){
    if( ut[i] ){
      if( is.na( x[i]) ){
        rowgroup <- df[ which(df$ID ==  rownames(x)[ row(x)[i] ]), "group" ]
        colgroup <- df[ which(df$ID ==  colnames(x)[ col(x)[i] ]), "group" ]
        z <- b[rowgroup, colgroup]
      } else z <- x[i]
    } else z <- NA
    return( z )
  } )

  vals <- as.numeric(unlist(out))

  y <- matrix( vals, length(ids), length(ids), dimnames=list(ids,ids) )
})
simmatrix <- y
saveRDS(simmatrix, file = file.path( wd$bin, "simmatrix.rds") )


# Perform hierarchical clustering -----------------------------------------

# For each species, perform normal clustering based on similarity matrixes. --
simmatrix <- readRDS(file.path( wd$bin, "simmatrix.rds"))
simmatrix_dist <- dist( 1 - simmatrix, diag = TRUE, upper = TRUE )
clustered_simmatrix <- hclust( simmatrix_dist, method = "ward.D2")
save(clustered_simmatrix, file = file.path(wd$bin, "clustered_simmatrices.Rdata") )


#  Prepare to cut trees into groups of individuals with common origins ----------------

# k-means cluster number
wss <- rep(NA, 15)
for (i in 1:15) wss[i] <- sum(kmeans(simmatrix_dist, i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# Plot trees.
plot(clustered_simmatrix, cex = 0.6)
rect.hclust(clustered_simmatrix, k = 6, border = 2:7)


# Estimate the latitude of highest probability of origin for each surface ------
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
maps_cropped_df    <- readRDS( file.path(wd$bin, "maps_df.rds") ) %>% 
  dplyr::filter(method == "raw")

ave_lat <- lapply(1:nrow(mydata_transformed), function(i){
  ex <- mydata_transformed[i, ]
  u1 <- na.omit( maps_cropped_df[ maps_cropped_df$ID == ex$ID , ] )
  lats <- base::sample(size = 1000, x = u1$y, replace = TRUE, prob = u1$value)
  mean(lats)
} ) %>%
  unlist

mydata_aveLat <- mydata_transformed %>% dplyr::mutate(ave_lat = ave_lat)

saveRDS(mydata_aveLat, file = file.path(wd$bin, "mydata_aveLat.rds"))

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
myClustDeets$ID <- mytree$labels

# Reorder cluster numbers with respect to descending latitude.
myClustDeets0 <- mydata_aveLat %>%
  right_join(myClustDeets, by = "ID")
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
  dplyr::select(ID, OriginCluster)

myClustDeets$originColor <- viridis::viridis(
  length(unique(myClustDeets$OriginCluster)))[myClustDeets$OriginCluster]

mydata_clustered <- myClustDeets

save(mydata_clustered, file = file.path(wd$bin, "mydata_clustered.Rdata"))

