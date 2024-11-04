library(terra)
library(tidyterra)
library(patchwork)

keratin_augsep    <- terra::rast(file.path(wd$bin, "keratin_augsep_isoscape.tif"))
keratin_augsep_se <- terra::rast( file.path(wd$bin, "iso_se_augsep.tif") )

iso <- c(keratin_augsep, keratin_augsep_se)
names(iso) <- c("iso", "iso_se")

# Plot isoscape
p1 <- ggplot() +
  geom_spatraster(iso[[1]], mapping = aes()) +
  scale_fill_viridis_c(
    expression(paste(delta^2~H[feather]," (", "\u2030", ", VSMOW-SLAP)") ),
    na.value = NA,
    option = "plasma", breaks = seq(-300,300,by=25)
    ) +
  scale_x_continuous(breaks = seq(-180,180,by=10)) +
  scale_y_continuous(breaks = seq(-180,180,by=10)) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

# Plot isoscape se
p2 <- ggplot() +
  geom_spatraster(iso[[2]], mapping = aes()) +
  scale_fill_viridis_c(
    expression(paste(delta^2~H[feather]," standard error (", "\u2030", ")") ),
    na.value = NA,
    option = "mako", breaks = c(0,5,10,15), limits = c(0,15)
  ) +
  scale_x_continuous(breaks = seq(-180,180,by=10)) +
  scale_y_continuous(breaks = seq(-180,180,by=10)) +
  theme_minimal()

# Combine
p_out <- p1/p2
ggsave(p_out, filename = file.path("figs", "SI_featherIsoscapes.png"), width = 6, height = 8)
