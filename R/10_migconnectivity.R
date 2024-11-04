library(MigConnectivity)
library(raster)
library(sf)
library(stringr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)

##### State map
state_prov <- ne_states(c("united states of america", "canada"))
state_prov <- st_as_sf(state_prov)

#### Locate isotope maps
f <- list.files("out/combined_maps")
f <- f[grepl("cont.tif", f)]
f <- f[!grepl("c_cont.tif", f)]

#### Create RasterStack from .tif files
IsoMaps <- stack(file.path("out/combined_maps", f))
names(IsoMaps) <- f
IsoMaps <- projectRaster(IsoMaps, crs = "+proj=longlat +datum=WGS84 +no_defs")


#### Extract band info for each map
meta <- data.frame(ID = word(f, start = 1,sep = "\\_"))

#### Match each map to origin lat/lon, and spp
xy <- read.csv("data/rosyFinch_dD_coords.csv")
meta <- left_join(meta, xy) %>% select(ID, Species, decimalLongitude, decimalLatitude)


##########################
### BLRF MC
##########################
blrf <- which(meta$Species == "BLRF")
blrfMeta <- meta[blrf,]
nblrf <- nrow(blrfMeta)
blrfMaps <- IsoMaps[[blrf]]

#### Save banding locations as sf object
originPoints <- sf::st_as_sf(data.frame(blrfMeta),
                             coords = c("decimalLongitude", "decimalLatitude"),
                             crs = sf::st_crs(state_prov))


isGL <-rep(F, nblrf); isRaster <- rep(T, nblrf)
isProb <- rep(F, nblrf); isTelemetry <- rep(F,nblrf)


originSite1 <- data.frame(lon = c(-112.25, -111),
                          lat = c(40, 41),
                          region = c("Powder", "Powder")) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(originPoints)) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

originSite2 <- data.frame(lon = c(-112.25, -111),
                          lat = c(41, 42),
                          region = c("Alta", "Alta")) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(originPoints)) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

originSite3 <- data.frame(lon = c( -110, -109),
                          lat = c(40.6, 41.25),
                          region = c("Dutch John", "Dutch John")) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(originPoints)) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

originSites <- dplyr::bind_rows(originSite1, originSite2, originSite3)

blrfTargetSites <- filter(state_prov, name %in% c("Utah", "Idaho", "Nevada","Montana", "Wyoming")) %>%
  dplyr::select(name, geometry)
blrfTargetSites <- st_cast(blrfTargetSites, to = "MULTIPOLYGON")

ggplot() +
  geom_sf(data = blrfTargetSites, color = "#2b2b2b", fill = "white", size=0.125) +
  geom_sf(data = originSites, color = "#2b2b2b", fill = "grey70", size=0.125) +
  geom_sf(data = originPoints) +
  coord_sf(crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"), datum = NA) +
  ggthemes::theme_map()

psi.blrf <- estTransition(isGL = isGL,
                      isRaster = isRaster,
                      isProb = isProb,
                      isTelemetry = isTelemetry,
                      targetSites = blrfTargetSites,
                      resampleProjection = sf::st_crs(OVENdata$targetPoints),
                      targetRaster = blrfMaps,
                      originSites = originSites,
                      originPoints = originPoints,
                      captured = rep("origin", nblrf),
                      targetNames = blrfTargetSites$name,
                      originNames = c("Power", "Alta", "Dutch John"),
                      verbose = 3,
                      nSamples = 100)
plot(psi.blrf)

originDist <- as.matrix(dist(st_distance(st_centroid(originSites))))
targetDist <- as.matrix(dist(st_distance(st_centroid(blrfTargetSites)), diag = T))

blrfMC <- estStrength(originDist = originDist,
             targetDist = targetDist,
             originRelAbund = c(0.08, 0.39, 0.53),
             psi = psi.blrf)


##########################
### GCRF MC
##########################
gcrf <- which(meta$Species == "GCRF")
gcrfMeta <- meta[gcrf,]
ngcrf <- nrow(gcrfMeta)
gcrfMaps <- IsoMaps[[gcrf]]

isGL <-rep(F, ngcrf); isRaster <- rep(T, ngcrf)
isProb <- rep(F, ngcrf); isTelemetry <- rep(F,ngcrf)

#### Save banding locations as sf object
originPoints <- sf::st_as_sf(data.frame(gcrfMeta),
                             coords = c("decimalLongitude", "decimalLatitude"),
                             crs = sf::st_crs(state_prov))


originSite1 <- data.frame(lon = c(-112.25, -111),
                          lat = c(40, 41),
                          region = c("1", "1")) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(originPoints)) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

originSite2 <- data.frame(lon = c(-112.25, -111),
                          lat = c(41, 42),
                          region = c("2", "2")) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = st_crs(originPoints)) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf()

originSites <- dplyr::bind_rows(originSite1, originSite2)

targetSites1 <- st_union(filter(state_prov, name %in% c("Washington", "Alberta", "British Columbia"))) %>%
  st_as_sf() %>%
  mutate(name = "Alberta Columbia")

st_geometry(targetSites1) <- "geometry"

targetSites2 <- st_union(filter(state_prov, name %in% c("California", "Oregon"))) %>%
  st_as_sf() %>%
  mutate(name = "USA")
st_geometry(targetSites2) <- "geometry"

targetSites3 <- filter(state_prov, name %in% c("Alaska", "Northwest Territories", "Yukon", "Washington", "Alberta", "British Columbia")) %>%
  dplyr::select(name, geometry)


gcrfTargetSites <- bind_rows(targetSites2, targetSites3)
gcrfTargetSites <- st_cast(gcrfTargetSites, to = "MULTIPOLYGON")

ggplot() +
  geom_sf(data = gcrfTargetSites, color = "#2b2b2b", fill = "white", size=0.125) +
  geom_sf(data = originSites, color = "#2b2b2b", fill = "grey70", size=0.125) +
  geom_sf(data = originPoints) +
  coord_sf(crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"), datum = NA) +
  ggthemes::theme_map()

psi.gcrf <- estTransition(isGL = isGL,
                      isRaster = isRaster,
                      isProb = isProb,
                      isTelemetry = isTelemetry,
                      targetSites = gcrfTargetSites,
                      # resampleProjection = sf::st_crs(OVENdata$targetPoints),
                      targetRaster = gcrfMaps,
                      originSites = originSites,
                      originPoints = originPoints,
                      captured = rep("origin", ngcrf),
                      targetNames = gcrfTargetSites$name,
                      originNames = c("Powder", "Alta"),
                      verbose = 3,
                      nSamples = 100)
plot(psi.gcrf)

originDist <- as.matrix(dist(st_distance(st_centroid(originSites))))
targetDist <- as.matrix(dist(st_distance(st_centroid(gcrfTargetSites)), diag = T))

gcrfMC <- estStrength(originDist = originDist,
            targetDist = targetDist,
            originRelAbund = c(0.47, 0.53),
            psi = psi.gcrf)

##### Difference in MC ----
diffMC(estimates = list(BLRF = blrfMC, GCRF = gcrfMC), nSamples = 1000)
