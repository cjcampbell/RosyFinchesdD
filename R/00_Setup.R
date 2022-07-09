
# Load packages.
# if(packageVersion("isocat") < "0.2.4.9000") warning("Install development version of isocat.")
# # devtools::install_github("cjcampbell/isocat")

# Load libraries required for analyses.
library(isocat)
library(dplyr)
library(sf)

# Load libraries relied on for a few analyses.
library(assignR)
library(chron)
library(geosphere)
library(ggmap)
library(ggplot2)
library(ggstar)
library(lubridate)
library(lwgeom)
library(measurements)
library(purrr)
library(readxl)
library(rgdal)
library(rmapshaper)
library(smatr)
library(stringr)
library(tidyr)

# Make an object to help navigate the subdirectories.
my_dir_path <-"/Users/cjcampbell/RosyFinchesdD"
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$data    <- file.path( my_dir_path, "data" )
wd$iucn    <- file.path( my_dir_path, "data", "iucn" )
wd$figs    <- file.path( my_dir_path, "figs" )
wd$out     <- file.path( my_dir_path, "out" )

# Check for presence of subdirectories. Create if needed.
invisible({
  lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) )
})


# Define extent of spatial analysis.
# Unit == meters
my_extent_aea <- raster::extent(
  -42e5, 18e5,
  -30e5, 32e5
)

# CRS for aea projection:
myCRS <- "+proj=aea +lon_0=-117 +lat_1=30 +lat_2=64 +lat_0=46 +datum=WGS84 +units=m +no_defs"
