# Load two csv files and combine -----

filenames <- c("RF-results 210119.xlsx", "Rosy-Finch_Utah_Results_211006.XLSX")
dat_loaded <- lapply( filenames, function(x) {
  whichfiles <- list.files( path = wd$data, pattern = x , full.names = TRUE)
  whichfile <- grep(whichfiles, pattern = "\\$", invert = T, value = T)
  thisdat   <- readxl::read_excel(whichfile)
  return(thisdat)
  }
  ) %>%
  dplyr::bind_rows()

# Load the files.
mydata <- dplyr::mutate(dat_loaded, ID = make.names(`Sample ID`))
