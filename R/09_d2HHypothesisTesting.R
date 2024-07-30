
source(file.path("R", "00_Setup.R"))
source(file.path("R", "01_loadIsotopeData.R"))

df <- mydata %>%
  dplyr::filter(sampleType == "Feather") %>%
  dplyr::mutate(
    subspecies = case_when(grepl("Hepburns", Notes) ~ "GCRF-Hep", TRUE ~ Species) ,
    yday = yday(date),
    location2 = case_when(grepl("Alta", Location) ~ "Alta", TRUE ~ Location)
    )

theme_set(theme_bw())
options(na.action = na.fail)


# Hypothesis testing ------------------------------------------------------

## Does d2H vary with respect to sex, month, etc?

ggplot(dplyr::filter(df, Species == "BLRF")) +
  geom_point(aes(y=d2H, x = date, color = location2, shape = Sex)) +
  scale_shape_manual(values = c(1, 2, 22)) +
  facet_wrap(~Species)

df1 <- dplyr::filter(df, Sex != "Unknown")

m1 <- lm(
  d2H ~ Location + Sex + yday,
  data = df1)
MuMIn::dredge(m1)
sjPlot::plot_model(m1)


df_female <- dplyr::filter(df2,  Species == "BLRF", Sex == "Female")
df_male <- dplyr::filter(df2,  Species == "BLRF", Sex == "Male")
t.test(df_female$d2H, df_male$d2H)


# Add in map results ------------------------------------------------------

if(!exists("distDir_ests")) distDir_ests <- fread( file.path(wd$bin, "distDir_ests.csv"))
if(!exists("minDistDeets")) minDistDeets <- fread(file.path(wd$bin, "minDistDeets"))
mydata_clustered <- readRDS(file.path(wd$bin, "mydata_clustered.rds"))

df2 <- distDir_ests %>%
  left_join(., mydata_clustered) %>%
  left_join(minDistDeets) %>%
  inner_join(df1)

dplyr::filter(df2, Species == "BLRF") %>%
  ggplot() +
  geom_point(aes(y=d2H, x = yday, color = dist_50, shape = location2)) +
  facet_wrap(~Species)
s
