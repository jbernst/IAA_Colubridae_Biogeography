# This script will extract occurrence data from the GBIF downloads.

# Set working directory
setwd("D:/Documents/Publications/Lycodon-Stegonotus/GBIF_obtained-data/")

# load packages
library(dplyr)

# list all of your zipped files in the current working directory
zip.files <- list.files(pattern = "\\.zip$")

# loop through each zipped file
for (zip.file in zip.files) {

  # extract the base name like "Genus_species"
  base.name <- sub("(_[^_]+)?\\.zip$", "", zip.file)
  
  # create temporary directory to place all the unzipped file contents into
  temp.dir <- tempfile()
  dir.create(temp.dir)
  
  # unzip the file into the temporary directory
  unzip(zip.file, exdir = temp.dir)
  
  # look for the "occurrence.txt" file inside our GBIF unzipped directory
  occurrence.file <- list.files(temp.dir, pattern = "occurrence.txt$", full.names = TRUE)

  # get the occurence records and merge the institutional code and catalog number, and also get the longitude and latitude, excluding NAs. Make this a new dataframe and write CSV
  if (length(occurrence.file) == 1) {
        occ <- read.delim(occurrence.file, stringsAsFactors = FALSE)
    if (all(c("institutionCode", "catalogNumber", "decimalLatitude", "decimalLongitude") %in% colnames(occ))) {
      df <- occ %>%
        mutate(
          ID = gsub(" ", "-", paste0(institutionCode, "_", catalogNumber)),
          lon = decimalLongitude,
          lat = decimalLatitude
        ) %>%
        filter(!is.na(lon),!is.na(lat)) %>%
        select(ID, lon, lat)
      

      output_filename <- paste0(base.name, "_GBIFrecords.csv")
      write.csv(df, file = output_filename, row.names = FALSE)
    } else {
      warning(paste("Missing required columns in", zip.file))
    }
  } else {
    warning(paste("No occurrence.txt found in", zip.file))
  }
  
  # zip the original directory back up
  files_to_zip <- list.files(temp.dir, full.names = TRUE, recursive = TRUE)
  zip(zip.file, files = files_to_zip, flags = "-r9Xj")
  
  # clean up temp directory
  unlink(temp.dir, recursive = TRUE)
}
