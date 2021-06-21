
# Find location of csv file we want to load.
which_file <- list.files( path = wd$data, pattern = "RF-results 210119.xlsx", full.names = TRUE)

# Load the file.
mydata <- readxl::read_excel( which_file ) %>%
  dplyr::mutate(ID = make.names(`Sample ID`))
