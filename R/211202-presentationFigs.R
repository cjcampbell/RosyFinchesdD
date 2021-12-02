library(viridisLite)
library(tidyverse)

noam <- readRDS(file.path(wd$bin, "NoAm.rds"))

iso2 <- raster::raster("/Users/cjcampbell/BigZaddyData/isoscapes/66100_caitjcampbell_Annual_NoAm_Map_H_1980_2010/66100_caitjcampbell_Annual_NoAm_Map_H_1980_2010/predkrig.tiff")

iso2_df <- projectRaster(iso2, crs= myCRS) %>%
  SDMetrics::surface2df()

if(!exists("my_isoscapes")) load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)

iso3 <- my_isoscapes$assignR_d2h_world$isoscape
iso3_df <- iso3 %>% SDMetrics::surface2df()

iso4 <- raster::raster("/Users/cjcampbell/BigZaddyData/isoscapes/83507_caitjcampbell_WNoAm_JulyAugust/83507_caitjcampbell_WNoAm_JulyAugust/predkrig.tiff") %>%
  projectRaster(crs = myCRS)
iso4_df <- iso4 %>% SDMetrics::surface2df()


ranges <- lapply( c("BLRF", "GCRF"), function(speciesCode){
  list.dirs(wd$iucn) %>%
    grep(pattern = speciesCode, value = T) %>%
    grep(pattern = "species_data", value = T) %>%
    rgdal::readOGR(dsn = ., layer = "data_0") %>%
    st_as_sf(crs = 4326) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
    st_make_valid() %>%
    dplyr::mutate(SEASONAL = factor(SEASONAL, levels = 1:3))
})

myCoords <- list(
  coord_sf(xlim = c(-4e6, 1e6), ylim = c(-2e6, 3.0e6))
)

# Plots.

library(ggstatsplot)
mydata %>%
  ggplot() +
  aes(d2H, y = Species, color = Species, group = interaction(Species, sampleType)) +
  geom_boxplot() +
  theme_pubclean()

boxplots <- mydata %>%
  mutate(Species_sample = paste(Species, sampleType, sep = "_")) %>%
  ggbetweenstats(
    data = ., y = d2H, x = Species_sample,
    ggplot.component	= list(
      scale_y_continuous( name = "d2H (VSMOW-SLAP)")
    )
  )

ggsave(boxplots, file = file.path(wd$figs, "boxplots.png"), width = 6, height = 4)
#






# Isoscapes
p_isoscape <- ggplot() +
  geom_tile(iso2_df, mapping = aes(x=x,y=y,fill=value, color = value)) +
  scale_fill_viridis_c(  "Annual precip. d2H\n(VSMOW-SLAP)" , option = "plasma" , limits = c(-227.004, -3.35529)) +
  scale_color_viridis_c( "Annual precip. d2H\n(VSMOW-SLAP)" , option = "plasma" , limits = c(-227.004, -3.35529)) +
  geom_sf(noam, mapping = aes(), fill = NA, size = 0.35) +
  myCoords +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  theme(
    plot.background = element_rect("white", color = "white"),
    panel.grid.major = element_line("grey80")
        ) +
  ggtitle("Annual isoscape")
ggsave(p_isoscape, file = file.path(wd$figs, "AnnualIsoscape.png"), width = 8, height = 4)

p_isoscape_gg <- ggplot() +
  geom_tile(iso3_df, mapping = aes(x=x,y=y,fill=value, color = value)) +
  scale_fill_viridis_c(  "Growing season precip. d2H\n(VSMOW-SLAP)" , option = "plasma",limits = c(-227.004, -3.35529)) +
  scale_color_viridis_c( "Growing season precip. d2H\n(VSMOW-SLAP)" , option = "plasma",limits = c(-227.004, -3.35529)) +
  geom_sf(noam, mapping = aes(), fill = NA, size = 0.35) +
  myCoords +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  theme(
    plot.background = element_rect("white", color = "white"),
    panel.grid.major = element_line("grey80")
  ) +
  ggtitle("Growing season isoscape")
ggsave(p_isoscape_gg, file = file.path(wd$figs, "GrowingSeasonIsoscape.png"), width = 8, height = 4)

p_isoscape_JA <- ggplot() +
  geom_tile(iso4_df, mapping = aes(x=x,y=y,fill=value, color = value)) +
  scale_fill_viridis_c(  "July/August precip. d2H\n(VSMOW-SLAP)" , option = "plasma",limits = c(-227.004, -3.35529)) +
  scale_color_viridis_c( "July/August precip. d2H\n(VSMOW-SLAP)" , option = "plasma",limits = c(-227.004, -3.35529)) +
  geom_sf(noam, mapping = aes(), fill = NA, size = 0.35) +
  myCoords +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  theme(
    plot.background = element_rect("white", color = "white"),
    panel.grid.major = element_line("grey80")
  ) +
  ggtitle("July/August isoscape")
ggsave(p_isoscape_JA, file = file.path(wd$figs, "JulyAugustIsoscape.png"), width = 8, height = 4)


p1 <- p_isoscape +
  ggtitle("BLRF range (IUCN)") +
  new_scale_color() +
  new_scale_fill() +
  geom_sf(ranges[[1]], mapping = aes(fill = SEASONAL, color = SEASONAL), alpha = 0.9, size = 1)+
  scale_fill_manual( "Seasonal range", values = c("#440154FF", "red", "blue"), labels = c("Year-round", "Breeding", "Nonbreeding"), drop = F) +
  scale_color_manual("Seasonal range", values = c("#440154FF", "red", "blue"), labels = c("Year-round", "Breeding", "Nonbreeding"), drop = F) +
  myCoords
ggsave(p1, file  = file.path(wd$figs, "BLRF_Range.png"), width = 8, height = 4)

p2 <- p_isoscape +
  ggtitle("GCRF range (IUCN)") +
  new_scale_color() +
  new_scale_fill() +
  geom_sf(ranges[[2]], mapping = aes(fill = SEASONAL, color = SEASONAL), alpha = 0.9, size = 1)+
  myCoords  +
  scale_fill_manual( "Seasonal range", values = c("#440154FF", "red", "blue"), labels = c("Year-round", "Breeding", "Nonbreeding")) +
  scale_color_manual("Seasonal range", values = c("#440154FF", "red", "blue"), labels = c("Year-round", "Breeding", "Nonbreeding"))
ggsave(p2, file  = file.path(wd$figs, "GCRF_Range.png"), width = 8, height = 4)


# Tiny b+w isoscape
colfun <- colorRampPalette(c("grey20", "grey70"))
mycols <- colfun(5)
p3 <- ggplot() +
  geom_tile(iso2_df, mapping = aes(x=x,y=y,fill=value, color = value)) +
  scale_color_stepsn(colors = mycols ) +
  scale_fill_stepsn( colors = mycols) +
  coord_sf() +
  theme_nothing()
p3
ggsave(p3, file = file.path(wd$figs, "exampleIsoscape.png"))
