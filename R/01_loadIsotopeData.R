# Load dD results  -----

filenames <- c("RF-results 210119.xlsx", "Rosy-Finch_Utah_Results_211006.XLSX")
dat_loaded <- lapply( filenames, function(x) {
  whichfiles <- list.files( path = wd$data, pattern = x , full.names = TRUE)
  whichfile <- grep(whichfiles, pattern = "\\$", invert = T, value = T)
  thisdat   <- readxl::read_excel(whichfile)
  return(thisdat)
  }
  ) %>%
  dplyr::bind_rows()


# Load metadata ----
dat_meta1 <-  list.files(
    path = wd$data, pattern = "RF-results 210119_w bird data.XLSX" ,
    full.names = TRUE
  ) %>%
    grep(., pattern = "\\$", invert = T, value = T)  %>%
    readxl::read_excel( sheet = "tbl_Temp") %>%
  dplyr::rename(
    date = `Date Banded`,
    OurLabID = `U of Maryland LabID`,
    `Sample ID` = `project Sample ID`,
    sampleType = `Type of Sample`
  )

dat_meta2 <-
  list.files(
    path = wd$data, pattern = "Rosy-Finch Keratin Analysis Submitted.xlsx" ,
    full.names = TRUE
  ) %>%
  grep(., pattern = "\\$", invert = T, value = T) %>%
  readxl::read_excel(
    sheet = "Combined",
    col_types = c("text", "date", "text", "text", "text", "text", "text",
                  "text", "text", "text", "text", "numeric", "text")
    ) %>%
  dplyr::rename(
    date = `Full Date Banded`,
    Location = `Location Banded/Sampled`,
    sampleType = `Sample Type`
  ) %>%
  # A few dates aren't read properly. Manual correction here:
  dplyr::mutate(
    date = case_when(
      is.na(date) & `Sample ID` %in% c("2441-20562c", "2441-20562f", "2441-20569c", "2441-20569f") ~ as.Date("02/05/2021", format = c("%m/%d/%Y")),
      TRUE ~ as.Date(date)
      )
    )

dat_meta <- full_join(dat_meta1, dat_meta2)


# Location data -------
locationData <- data.frame(
  Location = c("Alta", "Alta Town Office", "Powder Mountain", "Dutch John"),
  decimalLatitude = c(rep(40.590428, 2), 41.370703, 40.929955),
  decimalLongitude = c(rep(-111.637230,2 ), -111.768111, -109.389239)
)


# Combine --------
mydata <- full_join(dat_loaded, dat_meta) %>%
  dplyr::mutate(ID = make.names(`Sample ID`)  ) %>%
  dplyr::left_join(.,locationData) %>%
  dplyr::filter(!is.na(d2H))

