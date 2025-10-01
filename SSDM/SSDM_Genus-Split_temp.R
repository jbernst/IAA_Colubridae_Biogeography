# Note: this assumes you already have species observsation & environmental data already in hand, and also that you
# already have ran code to obtain rarefied datasets and your environmental variables ready.

# NOTE: This code will be updated upon publication of the manuscript.

##############################################
### Part 1: load packages
### Packages for clim data wranging and ENM
##############################################

setwd("D:/Documents/Publications/Lycodon-Stegonotus/SSDM/")

# load packages # list is extensive, may not need all
library(sf) # gis objects management
library(sp) # mapping
library(raster) # gis, raster management
library(terra) # gis, raster management (replaces package: raster)
library(maptools) # mapping
library(dismo) # sp data mining, modeling, gis
library(ecospat) # niche decomposition analysis
library(SSDM) # stacked niche models
library(heatmaply) # for correlation AnD heatmapping
library(tidyverse) # data wrangling
library(maps) # mapping
library(rJava) # java controls
#library(ENMeval) # enm tuning - I don't use it anymore
library(SDMtune) # enm tuning
#library(rasterVis) # raster management
#library(RColorBrewer) # color palettes
#library(phyloclim) # enm
#library(rgeos)
library(mapdata) # maps
library(gridExtra) # figure arrange
options(java.parameters = "-Xmx48g") # turn this on to give the maximum memory for Maxent niche modeling

# this assumed you have an external harddrive as an E: drive with a directory called R_temp in it. I did this 
# due to memory issues with the 30s resolution climate data, but this is not necessary (usually) with 2.5 m resolution.
# Though, the temp files R creates during the stacked models can be quite large. 
rasterOptions(tmpdir = "E:/R_temp")

##############################################
### Part 2: Climate data wrangler
##############################################

## 1. ## import climate data # CHANGE TO 2.5 MIN
raster_files <- list.files("D:/Documents/R/wc2.1_2.5m_bio/", full.names = T, pattern = ".tif")
  # if using WorldClim, remove variables with known discontinuities 
raster_files <- raster_files[-c(10,11)] # this removes bio_18 & bio_19 from the list above

# create a raster stack using the list you created
predictors <- stack(raster_files)
  # verify full clim data
plot(predictors$bio_1, main = "WC bio 1", xlab = "longitude", ylab = "latitude")

# get shapefiles (use this as a cookie cutter to trim environmental variables).
# if you ened to merge shape files, can use the package sf. Can also use terra (new version of the package raster) for
# typical raster management things.

MyCountryShape.shp <- st_read("D:/Documents/Mapping/DIVA-GIS_S-E-SEA-AA/DIVA-GIS_S-E-SEA-AA.shp")
  # verify 
plot(MyCountryShape.shp$geometry, main = "South Asia and the IAA") # must specify geometry for plotting only the polygon
  ## NOTE: sometimes shapefiles takes a ton of time to load depending on resolution and extent

# crop climate data to the calibration and/or projection extents
  # crop predictors to extent using shp file. This has two steps: crop and mask. Both steps are needed
Env <- crop(predictors, MyCountryShape.shp)
Env <- mask(Env, MyCountryShape.shp)
  # verify cropped data. it'll always look as if cropped with bounding box from shapefile
plot(Env$bio_1, main = "WC bio 1", xlab = "longitude", ylab = "latitude")

## Note:
## Don't forget to always verify each climate variable and look for climate discontinuities (abrupt zones of change)
  ## these discontinuities are typically shown in bio8, bio9, bio18, bio19 and cause hard cutoffs in distributions
  ## see Booth (2022) Austral Ecology 47, 1506-1514.

##############################################
### Part 3: Species data wrangler
##############################################

# Bring in rarefied dataset (see SSDM_All-Species_2025.05.29_Full.R on how to obtain these.
# observations should have 3 headers: species, lon, lat (x and y refer to lon and lat).

speciesNew_dir <- "D:/Documents/Publications/Lycodon-Stegonotus/SSDM/georeferenced_GBIF-records_GreaterThan10/rareFinal_datasets/"
spp.files <- list.files(speciesNew_dir, pattern = "\\.csv$", full.names = TRUE)

# you can assign each to an object
for (i in seq_along(spp.files)) {
  sp.rare.df <- read.csv(spp.files[i], header = TRUE)
  assign(paste0("species",i), sp.rare.df)
}

##############################################
### Part 4: Build SSDMs - split by genus
##############################################
# Check which genera are in which objects (note, R reads in the order of 1, 10, 11, 12, 13, so although the datasets are in generic order,
# the objects in R are not).

# head(species13) # LYCODON
# head(species14) # LYCODON
# head(species15) # LYCODON
# head(species16) # LYCODON
# head(species17) # LYCODON
# head(species18) # LYCODON

#=========#
# LYCODON #
#=========#

## 1. Import species dataframes

# use mget to retrieve multiple objects/files of our species dataframes (speceis1, species 12, etc...) and then combine them
# with rbind. 
Occ.lyco <- do.call(rbind, mget(paste0("species", 13:18)))  

# only take colums 1 to 3 (species, lon, lat)
Occ.lyco <- Occ.lyco[ , 1:3]

# confirm that only the genus and species we want is in the dataframe
unique(Occ.lyco$species)

# save the full community data set
write.csv(Occ.lyco, "D:/Documents/Publications/Lycodon-Stegonotus/SSDM/Occ_lyco.csv", row.names = F)

## 2. build SSDM using the respective genus

### A. stacked model using three algorithms MAXENT, GLM, and SVM ###
lyco.ssdm.mx.glm.svm <- stack_modelling(c('MAXENT', 'GLM', "SVM"), Occ.lyco, Env, rep = 1, ensemble.thresh = 0,
                                   Xcol = 'lon', Ycol = 'lat', Spcol = 'species', method = "pSSDM", verbose = T)

setwd("D:/Documents/Publications/Lycodon-Stegonotus/SSDM/")

pdf("Lyco_spp_ssdm.mx.glm.svm.pdf", lyco.ssdm.mx.glm.svm)
plot(lyco.ssdm.mx.glm.svm@diversity.map, main = 'Stacked Model for Lycodon')
dev.off()

## 3. obtain model evaluations and statistics, and print them to files

# get general evaluation of stacked model
sink("GeneralEval_Lycodon.txt")
knitr::kable(lyco.ssdm.mx.glm.svm@evaluation)
sink()

# get model statistics -- general evaluation
knitr::kable(lyco.ssdm.mx.glm.svm@evaluation)

sink("GeneralEval_Lycodon.txt")
print(knitr::kable(lyco.ssdm.mx.glm.svm@evaluation))
sink()

# get model statistics -- variable importance
knitr::kable(lyco.ssdm.mx.glm.svm@variable.importance)

sink("VariableImportance_Lycodon.txt")
print(knitr::kable(lyco.ssdm.mx.glm.svm@variable.importance))
sink()

# get model statistics -- algorithm evaluation
knitr::kable(lyco.ssdm.mx.glm.svm@algorithm.evaluation)

sink("AlgorithmEvaluation_Lycodon.txt")
print(knitr::kable(lyco.ssdm.mx.glm.svm@algorithm.evaluation))
sink()

# get model statistics -- algorithm correlation
knitr::kable(lyco.ssdm.mx.glm.svm@algorithm.correlation)

sink("AlgorithmCorrelation_Lycodon.txt")
print(knitr::kable(lyco.ssdm.mx.glm.svm@algorithm.correlation))
sink()

## 4. write out rasters of the SSDM

# richness map
writeRaster(lyco.ssdm.mx.glm.svm@diversity.map[[1]], filename = "D:/Documents/Publications/Lycodon-Stegonotus/SSDM/Lycodon_RichnessRaster.tif")

# endemism map
writeRaster(lyco.ssdm.mx.glm.svm@endemism.map[[1]], filename = "D:/Documents/Publications/Lycodon-Stegonotus/SSDM/Lycodon_EndemismRaster.tif")
