# Load  ----------------------------------------------------------
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopedata.R"))

library(assignR)
if(is.null(wd$isoscapes)) warning("Assign wd$isoscapes as the location of the
                                    isoscapes on the local computer.")

# Load global isoscapes. -------------------------------------------------------
iso_GS_WGS84     <- raster::raster( file.path(wd$isoscapes, "GlobalPrecipGS", "d2h_GS.tif"))

# Load assignR known-origin dataset. -------------------------------------------
assignR_Hobson <- assignR::knownOrig$samples %>%
  dplyr::filter(Dataset_ID == 2) %>%
  dplyr::inner_join(., as.data.frame(assignR::knownOrig$sites), by = "Site_ID") %>%
  dplyr::inner_join(., as.data.frame(assignR::knownOrig$sources), by = "Dataset_ID")

# Use refTrans function to translate to VSMOW reference scale.
assignR_Hobson_VSMOW <- pbmcapply::pbmclapply( 1:nrow(assignR_Hobson), mc.cores = 3, function(i) {
  df <- assignR_Hobson[i, ]
  s <- data.frame("d2H" = df$d2H, "d2H.sd" = df$d2H.sd, "d2H_cal" = df$H_cal)
  df1 <- assignR::refTrans(s)
  out <- data.frame(df, df1$data)
  return(out)
}) %>% bind_rows()

# Plot changes
assignR_Hobson_VSMOW %>%
  ggplot() +
  geom_point(aes(x=d2H, y = d2H.1)) +
  geom_abline(slope=1,intercept = 0, color = "blue")

# Find dPrecip values at sample sites ------------------------------------------
assignR_Hobson_sf <- assignR_Hobson_VSMOW %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

assignR_Hobson2 <- dplyr::mutate(
  assignR_Hobson_VSMOW,
  iso_GS_WGS84 = raster::extract(iso_GS_WGS84, assignR_Hobson_sf)
)

# Exploratory plot.
assignR_Hobson2 %>%
  ggplot() +
  geom_point(aes(x=iso_GS_WGS84, y = d2H.1))


# Bring in guild predictors from Hobson 2012 SI --------------------------------

Hobson2012_SI1 <- read.csv(file.path(wd$data, "Hobson2012_SI1.csv") ) %>%
  dplyr::filter(!is.na(Common.Name), Common.Name != "") %>%
  dplyr::select(-n)

# Find names without a perfect match.
unique(assignR_Hobson2$Taxon)[!unique(assignR_Hobson2$Taxon) %in% Hobson2012_SI1$Scientific.Name]
unique(Hobson2012_SI1$Scientific.Name)[!unique(Hobson2012_SI1$Scientific.Name) %in% assignR_Hobson2$Taxon]

key <- data.frame(
  assignR_name = c(
    "Setophaga coronata auduboni", # Subspecies prevents perfect match
    "Junco hyemalis oregonus",     # Subspecies prevents perfect match
    "Vermivora cyanoptera",        # Synonym 'Vermivora pinus'
    "Poecile carolinensis"         # SI has a typo
  ),
  Hobson2012_name = c(
    "Setophaga coronata",
    "Junco hyemalis",
    "Vermivora pinus",
    "Poecile carlinensis"
  )
)
# Loggerheaded shrike in assignR ("Lanius ludovicianus") were not in the original
# Hobson 2012 data, though they were pooled in with Hobson-sourced data
# when assembling the assignR dataset.
#
# There are 30 tree swallows ("Tachycineta bicolor") that were included in the
# Hobson 2012 analyses that were not in the assignR dataset.


assignR_Hobson3 <- assignR_Hobson2 %>%
  dplyr::filter(Taxon != "Lanius ludovicianus") %>%
  left_join(., key, by = c("Taxon" = "assignR_name")) %>%
  dplyr::mutate(Hobson2012_name = case_when(
    is.na(Hobson2012_name) ~ Taxon,
    TRUE ~ Hobson2012_name
  )) %>%
  left_join(., Hobson2012_SI1, by = c("Hobson2012_name" = "Scientific.Name")) %>%
  dplyr::mutate(Foraging.Substrate = factor(Foraging.Substrate, levels = c("Non-Ground", "Ground")))
saveRDS(assignR_Hobson3, file = file.path(wd$bin, "assignR_Hobson3.rds"))

# Fit mixed models as in Hobson 2012. ------------------------------------------

m1 <- lm(d2H.1~iso_GS_WGS84 + Foraging.Substrate + Migratory + Foraging.Substrate*Migratory,
     data = assignR_Hobson3 )
summary(m1)
AIC(m1)
(a <- sjPlot::plot_model(m1, type = "pred", terms = c("iso_GS_WGS84", "Foraging.Substrate", "Migratory")))

# ForagingSubstrateGround == 1
# MigratoryShort distance == 1
# Foraging.SubstrateGround:MigratoryShort distance == 1
intercept <- coef(m1)[
  names(coef(m1)) %in% c(
    "(Intercept)",
    "Foraging.SubstrateGround",
    "MigratoryShort distance",
    "Foraging.SubstrateGround:MigratoryShort distance")
] %>%
  sum()
slope <- coef(m1)[2]
groundShort <- assignR_Hobson3 %>%
  dplyr::filter(Foraging.Substrate == "Ground", Migratory == "Short distance")
groundShort %>%
  ggplot() +
  aes(x=iso_GS_WGS84, y=d2H.1, color=Migratory, shape = Foraging.Substrate) +
  geom_point() +
  geom_abline(slope = slope, intercept = intercept)

params <- list(slope=slope, intercept=intercept)
saveRDS(params, file = file.path(wd$bin, "transferFunctionParams.rds"))



# Estimate sd  of resids------------------------------------------------------------
calculateResidual <- function(x,y, slope = slope, intercept = intercept) {
  yhat = (x*slope)+intercept
  resid = y-yhat
  return(resid)
}
resids <- mapply(FUN = calculateResidual, x = groundShort$iso_GS, y = groundShort$d2H.1,
       MoreArgs = list( slope = slope, intercept = intercept)
)
hist(resids)
sd_resids <- sd(resids)

saveRDS(sd_resids, file = file.path(wd$bin, "transferFunctionResids.rds"))


# Transform isoscape ------------------------------------------------------
precipToKeratin <- function(x) {
  stopifnot(exists("params"))
  y <- ( x*params$slope) + params$intercept
  return(y)
}

## Feather isoscape ----
# August-September surfaces
iso_augsep    <- raster::raster( file.path(wd$bin, "iso_augsep.tif") )
keratin_augsep <- raster::calc(iso_augsep, fun =  precipToKeratin)
writeRaster(keratin_augsep, filename = file.path(wd$bin, "keratin_augsep_isoscape.tif"), overwrite = T)

### Sep-Dec ----
iso_sepdec    <- raster::raster( file.path(wd$bin, "iso_sepdec.tif") )
keratin_sepdec <- raster::calc(iso_sepdec, fun =  precipToKeratin)
writeRaster(keratin_sepdec, filename = file.path(wd$bin, "keratin_sepdec_isoscape.tif"), overwrite = T)

### Oct-Jan ----
iso_octjan    <- raster::raster( file.path(wd$bin, "iso_octjan.tif") )
keratin_octjan <- raster::calc(iso_octjan, fun =  precipToKeratin)
writeRaster(keratin_octjan, filename = file.path(wd$bin, "keratin_octjan_isoscape.tif"), overwrite = T)

### Nov-Feb ----
iso_novfeb    <- raster::raster( file.path(wd$bin, "iso_novfeb.tif") )
keratin_novfeb <- raster::calc(iso_novfeb, fun =  precipToKeratin)
writeRaster(keratin_novfeb, filename = file.path(wd$bin, "keratin_novfeb_isoscape.tif"), overwrite = T)
