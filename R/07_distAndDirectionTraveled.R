

# distDirSummaryPlot ------------------------------------------------------


if(!exists("distDir_ests")) distDir_ests <- fread( file.path(wd$bin, "distDir_ests.csv"))
if(!exists("minDistDeets")) minDistDeets <- fread(file.path(wd$bin, "minDistDeets"))

p_dist_BLRF <- distDir_ests %>%
  left_join(., mydata_clustered) %>%
  left_join(minDistDeets) %>%
  dplyr::filter(Species == "BLRF") %>%
  dplyr::mutate( subspecies = case_when(grepl("Hepburns", Notes) ~ "GCRF-Hep", TRUE ~ Species) ) %>%
  dplyr::filter(sampleType == "Feather") %>%
  group_by(Species, sampleType) %>%
  arrange(desc(d2H)) %>%
  dplyr::mutate(rn = row_number()) %>%
  ggplot() +
  geom_rect(aes(xmin = dist_25, xmax = dist_75, ymin = rn - 0.3 , ymax = rn + 0.3 , color = subspecies), fill = NA) +
  geom_segment(aes(y = rn - 0.3 , yend = rn + 0.3 , x = dist_50, xend = dist_50, color = subspecies), linewidth = 2) +
  geom_segment(aes(y = rn, yend = rn, x = dist_05, xend =dist_25 ,color = subspecies)) +
  geom_segment(aes(y = rn, yend = rn, x = dist_75, xend =dist_95 ,color = subspecies)) +
  scale_x_continuous("Distance from origin to sample site (km)", limits = c(0,NA)) +
  scale_color_manual(values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  scale_fill_manual( values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  ylab(NULL) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x.top = element_line(),
    strip.background = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    legend.background = element_rect(fill = NA)
  )


p_dist_GCRF <- distDir_ests %>%
  left_join(., mydata_clustered) %>%
  left_join(minDistDeets) %>%
  dplyr::filter(Species == "GCRF") %>%
  dplyr::mutate( subspecies = case_when(grepl("Hepburns", Notes) ~ "GCRF-Hep", TRUE ~ Species) ) %>%
  dplyr::filter(sampleType == "Feather") %>%
  group_by(Species, sampleType) %>%
  arrange(desc(d2H)) %>%
  dplyr::mutate(rn = row_number()) %>%
  ggplot() +
  geom_rect(aes(xmin = dist_25, xmax = dist_75, ymin = rn - 0.3 , ymax = rn + 0.3 , color = subspecies), fill = NA) +
  geom_segment(aes(y = rn - 0.3 , yend = rn + 0.3 , x = dist_50, xend = dist_50, color = subspecies), linewidth = 2) +
  geom_segment(aes(y = rn, yend = rn, x = dist_05, xend =dist_25 ,color = subspecies)) +
  geom_segment(aes(y = rn, yend = rn, x = dist_75, xend =dist_95 ,color = subspecies)) +
  scale_x_continuous("Distance from origin to sample site (km)", limits = c(0,NA)) +
  scale_color_manual(values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  scale_fill_manual( values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  ylab(NULL) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x.top = element_line(),
    strip.background = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    legend.background = element_rect(fill = NA)
  )


p_direction_BLRF <- minDistDeets %>%
  dplyr::mutate(
    subspecies = case_when(grepl("Hepburns", Notes) ~ "GCRF-Hep", TRUE ~ Species)
  ) %>%
  dplyr::filter(Species == "BLRF", sampleType == "Feather") %>%
  ggplot() +
  geom_histogram(aes(theta_from_site, fill = subspecies), binwidth = 45/2) +
  coord_polar(theta = "x", start = pi) +
  geom_star(aes(x=0,y=0), size = 2, fill = "white") +
  scale_x_continuous(
    name = NULL,
    breaks = seq(-180,179,90),
    limits = c(-180,180),
    expand = c(0,0),
    labels = c("S", "W", "N", "E")
  ) +
  scale_fill_manual( values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x.top = element_line(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = NA),
    legend.position = "none"
  )

p_direction_GCRF <- minDistDeets %>%
  dplyr::mutate(
    subspecies = case_when(grepl("Hepburns", Notes) ~ "GCRF-Hep", TRUE ~ Species)
  ) %>%
  dplyr::filter(Species == "GCRF", sampleType == "Feather") %>%
  ggplot() +
  geom_histogram(aes(theta_from_site, fill = subspecies), binwidth = 45/2) +
  coord_polar(theta = "x", start = pi) +
  scale_x_continuous(
    name = NULL,
    breaks = seq(-180,179,90),
    limits = c(-180,180),
    expand = c(0,0),
    labels = c("S", "W", "N", "E")
  ) +
  geom_star(aes(x=0,y=0), size = 2, fill = "white") +
  scale_fill_manual( values = c("BLRF" = "grey20", "GCRF" = "grey30", "GCRF-Hep" = "darkgoldenrod"))+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x.top = element_line(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = NA),
    legend.position = c(0.10,0.1),
    legend.title = element_blank()
  )

ggarrange(p_dist_BLRF,p_direction_BLRF, p_dist_GCRF,p_direction_GCRF, ncol = 2, nrow = 2, labels = LETTERS) %>%
  ggsave(., filename = file.path(wd$figs, "distDIR-feathers_species.png"), width = 6, height = 6)

# Summarize distances traveled. ---------------
library(ggridges)
if(!file.exists(file.path(wd$bin, "distDir_ests.csv"))) {
  set.seed(42)
  distDir_ests <- mydat_distDir0 %>%
    # dplyr::filter(
    #   Species == "BLRF",
    #   sampleType == "Feather"
    # ) %>%
    left_join(., mydata_clustered) %>%
    dplyr::filter(
      # Species == "BLRF",
      # sampleType == "Feather",
      #ID %in% c("X2441.20646f", "X2441.20654f", "X2441.20656f", "X2441.20663f", "X2441.20672f", "X2441.20673f"),
      !is.na(combo1)
    ) %>%
    group_by(ID) %>%
    sample_n(size = 10000, weight = combo1, replace = T) %>%
    dplyr::summarize(
      dist_05 = quantile(dist_km, 0.05),
      dist_25 = quantile(dist_km, 0.25),
      dist_50 = quantile(dist_km, 0.50),
      dist_75 = quantile(dist_km, 0.75),
      dist_95 = quantile(dist_km, 0.95),
      dir_mean = mean(theta_from_site)
    )
  fwrite(distDir_ests, file = file.path(wd$bin, "distDir_ests.csv"), row.names = F)
}
if(!exists("distDir_ests")) distDir_ests <- fread( file.path(wd$bin, "distDir_ests.csv"))
