source("~/RosyFinchesdD/.Rprofile")

# Setup -------------------------------------------------------------------
library(patchwork)
library(ggpubr)

if(!exists("mydat_distDir0")) {  mydat_distDir0 <- fread(file.path(wd$bin, "mydat_distDir0")) }
if(!exists("minDistDeets")) {  minDistDeets <- fread(file.path(wd$bin, "minDistDeets")) }
states <- readRDS(file.path(wd$bin, "states.rds"))
if(!exists("elev")) {
  elev0 <- list.files(wd$bioclim, pattern = "elev.tif", full.names = T, recursive = T) %>%
    terra::rast() %>%
    terra::project(myCRS) %>%
    terra::crop(my_extent_aea) %>%
    raster::raster()
  elev <- elev0 %>%
    SDMetrics::surface2df()
}

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
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 8),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    title = element_text(size = 13),
    legend.direction = "horizontal",
    legend.title.align = 0.5
  )
)



# Plot maps ---------------------------------------------------------------

myIDs <- minDistDeets %>%
  data.frame %>%
  dplyr::filter(sampleType == "Feather") %>%
  group_by(Species) %>%
  dplyr::mutate(diff = abs(d2H - median(d2H))) %>%
  arrange(diff) %>%
  slice_sample(n=1) %>%
  ungroup %>%
  dplyr::select(ID) %>%
  unlist

myplots <- list()
for(myID in myIDs) {

  mySpecies <- unlist(minDistDeets[minDistDeets$ID == myID, "Species"])

  if(mySpecies == "BLRF") {
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

  mdf <- dplyr::filter(mydat_distDir0, ID == myID)
  p1_d2H <- ggplot() +
    geom_tile(
      data = mdf,
      aes(x=origin_x_m,y=origin_y_m,color=value, fill=value)
    ) +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    geom_star( data = dplyr::filter(minDistDeets, ID == myID), aes(x=sample_x_m, y= sample_y_m), size = 2, fill = "white") +
    moreDeets +
    scale_fill_viridis_c(  "Probability of origin", option = "mako", na.value = NA, breaks = range(mdf$value, na.rm = T), labels = c("Low", "High")) +
    scale_color_viridis_c( "Probability of origin", option = "mako", na.value = NA, breaks = range(mdf$value, na.rm = T), labels = c("Low", "High"))

  p1_abund <- ggplot() +
    geom_tile(
      data = mdf,
      aes(x=origin_x_m,y=origin_y_m,color=abund, fill=abund)
    ) +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    geom_star( data = dplyr::filter(minDistDeets, ID == myID), aes(x=sample_x_m, y= sample_y_m), size = 2, fill = "white") +
    moreDeets +
    scale_fill_viridis_c(  "Predicted abundance", option = "plasma", na.value = NA, breaks = range(mdf$abund, na.rm = T), labels = c("Low", "High")) +
    scale_color_viridis_c( "Predicted abundance", option = "plasma", na.value = NA, breaks = range(mdf$abund, na.rm = T), labels = c("Low", "High"))

  p1_combo <- ggplot() +
    geom_tile(
      data = mdf,
      aes(x=origin_x_m,y=origin_y_m,color=combo2_OR, fill=combo2_OR)
    ) +
    geom_sf(antiPoly, mapping = aes(), fill = "grey40", color = NA) +
    plotDeets +
    geom_star( data = dplyr::filter(minDistDeets, ID == myID), aes(x=sample_x_m, y= sample_y_m), size = 2, fill = "white") +
    moreDeets +
    scale_fill_viridis_c(  "Odds of origin", option = "turbo", na.value = NA, breaks = seq(0,1,by = .25), limits = c(0,1)) +
    scale_color_viridis_c( "Odds of origin", option = "turbo", na.value = NA, breaks = seq(0,1,by = .25), limits = c(0,1))


  myplots[[length(myplots) + 1]] <- p1_d2H
  myplots[[length(myplots) + 1]] <- p1_abund
  myplots[[length(myplots) + 1]] <- p1_combo
}

#
# library(gridExtra)
# pp <- arrangeGrob(
#   grobs = myplots,
#   layout_matrix = matrix(nrow = 2, ncol = 4, byrow = T, data = c(1,2,4,5,3,3,6,6)),
#   heights = c(1,3), labels = LETTERS[1:6])
# ggsave(pp, filename = file.path(wd$figs, paste0("median_combo.png")), width = 6, height = 6)
#
#
#
# library(patchwork)
#
# (myplots[[1]] + myplots[[4]] + plot_layout(guides = "collect"))
#
# library(ggpubr)
# ppp <- ggarrange(
#   ggarrange(
#     myplots[[1]] + theme(legend.position = "bottom"),
#     myplots[[2]] + theme(legend.position = "bottom"), ncol = 2,
#     labels = LETTERS[1:2]
#   ),
#
#   ggarrange(
#     myplots[[4]] + theme(legend.position = "bottom"),
#     myplots[[5]] + theme(legend.position = "bottom"), ncol = 2,
#     labels = LETTERS[4:5]
#   ),
#
#   ggarrange(myplots[[3]] + theme(legend.position = "bottom"), labels = LETTERS[3]),
#   ggarrange(myplots[[6]] + theme(legend.position = "bottom"), labels = LETTERS[6]),
#   ncol = 2, nrow = 2, heights = c(1,2)
# )
# ggsave(ppp, filename = file.path(wd$figs, paste0("median_combo2.png")), width = 4, height = 8)


l1 <- get_legend(myplots[[1]] + theme(legend.text = element_text(size = 5), legend.title = element_text(size = 5), legend.key.width = unit(.33, "cm"), legend.key.height = unit(.20, "cm")))
l2 <- get_legend(myplots[[2]] + theme(legend.text = element_text(size = 5), legend.title = element_text(size = 5), legend.key.width = unit(.33, "cm"), legend.key.height = unit(.20, "cm")))
l3 <- get_legend(myplots[[3]] + theme(legend.text = element_text(size = 5), legend.title = element_text(size = 5), legend.key.width = unit(.50, "cm"), legend.key.height = unit(.25, "cm")))


p1 <- cowplot::ggdraw() +
  cowplot::draw_plot(myplots[[1]] + theme(legend.position = "none"), x = 0.0, y = 1, width = 0.5, height = 1/3, hjust = 0, vjust = 1) +
  cowplot::draw_plot(myplots[[2]] + theme(legend.position = "none"), x = 1.0, y = 1, width = 0.5, height = 1/3, hjust = 1, vjust = 1) +
  cowplot::draw_plot(myplots[[3]] + theme(legend.position = "none"), x = 0.0, y = 0, width = 1.0, height = 2/3, hjust = 0, vjust = 0) +
  cowplot::draw_plot(l1, x = 0.02, y = .69, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.75) +
  cowplot::draw_plot(l2, x = 0.52, y = .69, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.75) +
  cowplot::draw_plot(l3, x = 0.10, y = 0.0, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.70) +
  cowplot::draw_label("A", x = 0.0, y = 1.0, hjust = 0, vjust = 1.0, size = 10) +
  cowplot::draw_label("B", x = 0.5, y = 1.0, hjust = 1, vjust = 1.0, size = 10) +
  cowplot::draw_label("C", x = 0.0, y = 2/3, hjust = 0, vjust = 0.5, size = 10)
p2 <- cowplot::ggdraw() +
  cowplot::draw_plot(myplots[[4]] + theme(legend.position = "none"), x = 0.0, y = 1, width = 0.5, height = 1/3, hjust = 0, vjust = 1) +
  cowplot::draw_plot(myplots[[5]] + theme(legend.position = "none"), x = 1.0, y = 1, width = 0.5, height = 1/3, hjust = 1, vjust = 1) +
  cowplot::draw_plot(myplots[[6]] + theme(legend.position = "none"), x = 0.0, y = 0, width = 1.0, height = 2/3, hjust = 0, vjust = 0) +
  cowplot::draw_plot(l1, x = 0.02, y = .69, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.75) +
  cowplot::draw_plot(l2, x = 0.52, y = .69, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.75) +
  cowplot::draw_plot(l3, x = 0.10, y = 0.0, width = 0.07, height = 0.2, hjust = -0.75, vjust = 0.70) +
  cowplot::draw_label("D", x = 0.0, y = 1.0, hjust = 0, vjust = 1.0, size = 10) +
  cowplot::draw_label("E", x = 0.5, y = 1.0, hjust = 1, vjust = 1.0, size = 10) +
  cowplot::draw_label("F", x = 0.0, y = 2/3, hjust = 0, vjust = 0.5, size = 10)

ggarrange(p1, p2, ncol = 2) %>%
  ggsave(filename = file.path(wd$figs, paste0("median_combo3.png")), width = 6, height = 4)

