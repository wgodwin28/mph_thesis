#Prep Mexico Shapefile
rm(list=ls())
#set up smart sensing of os
if(Sys.info()[1]=='Windows'){
  j = "J:/"
} else{
  j = '/home/j/'
}

# load packages, install if missing  
pack_lib = '/snfs2/HOME/wgodwin/R'
.libPaths(pack_lib)
pacman::p_load(data.table, fst, ggplot2, parallel, magrittr, maptools, raster, rgdal, rgeos, sp, splines, stringr, RMySQL, snow, ncdf4)

#Set toggles
iso <- "mex"
start <- 1987
end <- 1988
pop.dir <- paste0(j, "WORK/11_geospatial/01_covariates/00_MBG_STANDARD/worldpop/total/1y/")
lights.dir <- paste0(j, "WORK/11_geospatial/01_covariates/00_MBG_STANDARD/dmspntl/mean/1y/")
shapefile.dir = paste0(j, "WORK/05_risk/risks/temperature/data/exp/shapes/", iso, "/")
out.dir <- paste0(j, "temp/wgodwin/thesis/shapes/")
dir.create(cov.dir)
source(paste0(j, "temp/central_comp/libraries/current/r/get_ids.R"))

########Read in shapefile and extract##########
borders <- readOGR(paste0(shapefile.dir, iso, "_admin2.shp"))

# all mexico municipalities in mexico city
mex.city.munis <- c(9002, 9014, 9003, 9004, 9015, 9005, 9006, 9007, 9008, 9016, 9009, 9011, 9012, 9017, 9013, 9010)

#Subset shapefile
borders2 <- borders[borders@data$adm2_code %in% mex.city.munis,]
writeOGR(borders2, out.dir, "mex_city", driver="ESRI Shapefile")
