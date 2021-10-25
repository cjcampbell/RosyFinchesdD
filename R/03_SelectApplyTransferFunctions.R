# Load  ----------------------------------------------------------

if(!exists("my_isoscapes")) load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopedata.R"))


# Apply transfer functions, define molt status  ---------------------------------

# FUNCTION currently from Hobson 2021
# THIS IS A PLACEHOLDER AND ONLY WORKS FOR FEATHERS
# d2Hf = 0.95*d2Hp - 27.09
# MADE UP SD

mydata_transformed <- mydata %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::mutate(
    dDprecip =  (d2H - (-27.09)) / 0.95,
    sdResid = 5
  )

saveRDS(mydata_transformed, file = file.path(wd$bin, "mydata_transformed.rds"))


# Pull out selected isoscape ----------------------------------------------

bestFitIso <- my_isoscapes[[sma_selected$isoscape]]
save(bestFitIso, file = file.path(wd$bin, "bestFitIso.rdata"))

