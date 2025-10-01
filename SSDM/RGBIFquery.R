# rGBIF: hardcore query
# source: https://github.com/ropensci/rgbif
# edited by: Justin M. Bernstein

library(rgbif)

#===LYCODON FASCIATUS===#

# search species by name
pfal_rgbif <- occ_search(scientificName = "Lycodon fasciatus")
pfal_rgbif1 <- occ_data(scientificName = "Lycodon fasciatus")

# using taxon key for searches
taxonKey <- name_backbone("Lycodon fasciatus")$usageKey
occ_search(taxonKey = taxonKey)
# download records for [Genus species] (use the taxonKey): this is recommended to get DOI but it is more involved

# a. set up user attributes and account with gbif
# usethis::edit_r_environ()
# follow this info to add to R environment: https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
# the R environment must contain the following info:
# GBIF_USER="user_name"
# GBIF_PWD="user_password"
# GBIF_EMAIL="user@email.com"
# just need to do this once and restart R to ave your R environment

# b. query and download
occ_download(pred("taxonKey", 5222796)) # use taxon key given by the output of the above code
## processing takes a while and results be available for download manually
# see the following to check the status (the download key value is given from the above occ_download output: 
occ_download_wait('0023673-250426092105405')

# retrieve download from:
d <- occ_download_get('0023673-250426092105405') %>%
  occ_download_import()

# <<gbif download metadata>>
# Status: SUCCEEDED
# DOI: 10.15468/dl.yqqznn
# Format: DWCA
# Download key: 0023673-250426092105405
# Created: 2025-05-09T18:11:49.920+00:00
# Modified: 2025-05-09T18:13:32.233+00:00
# Download link: https://api.gbif.org/v1/occurrence/download/request/0023673-250426092105405.zip
# Total records: 204

# Repeat this with other taxa!
