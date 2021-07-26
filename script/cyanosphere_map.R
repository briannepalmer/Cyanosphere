# make map of cyanosphere points

# load libraries 
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(tidyverse)

cyanosphere <- read.csv("data/autometa_assembly_cyanosphere.csv")

world <- map_data("world") 
gg1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "gray") + 
  coord_fixed(1.3) + theme_bw()

gg1 + 
  geom_point(data = cyanosphere, aes(x = Long, y = Lat, color = Substrate), size = 1.5) + geom_jitter()

# add climate data
library(raster)

# run to get climate data
#climate <- getData('worldclim', var='bio', res=2.5)

#### Bio 1 ####
latlong <- cyanosphere %>% group_by(Long, Lat, culture_ID) %>% summarize()

point <- SpatialPoints(coords = cbind(latlong$Long,latlong$Lat))
crs <- crs(meantemp)

points_spdf <- SpatialPointsDataFrame(latlong[,1:2], proj4string= crs, latlong)

# bio1.df <- raster::extract(climate$bio1,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

# mean temp is C degrees times ten 

bio1.df$bio1 <- bio1.df$bio1/10

plot(climate$bio1, main = "Annual Mean Temperature") 
points(x = latlong$Long, y = latlong$Lat)


#### Bio 2 ####
# bio2.df <- raster::extract(climate$bio2,             # raster layer
#                                points_spdf,   # SPDF with centroids for buffer
#                                buffer = 1,     # buffer size, units depend on CRS
#                                fun=mean,         # what to value to extract
#                                df=TRUE)         # return a dataframe? 


plot(climate$bio2, main = "Mean Diurnal Range (Mean of monthly (max temp - min temp))") 
points(x = latlong$Long, y = latlong$Lat)


#### Bio 3 ####
# bio3.df <- raster::extract(climate$bio3,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 
# 
# plot(climate$bio3, main = "Isothermality (BIO2/BIO7) (×100)") 
# points(x = latlong$Long, y = latlong$Lat)

#### Bio 4 ####
# bio4.df <- raster::extract(climate$bio4,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio4, main = "Temperature Seasonality (standard deviation ×100)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 5 ####
# bio5.df <- raster::extract(climate$bio5,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio5, main = "Max Temperature of Warmest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 6 ####
# bio6.df <- raster::extract(climate$bio6,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio6, main = "Min Temperature of Coldest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 7 ####
# bio7.df <- raster::extract(climate$bio7,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio7, main = "Temperature Annual Range (BIO5-BIO6)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 8 ####
# bio8.df <- raster::extract(climate$bio8,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio8, main = "Mean Temperature of Wettest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 9 ####
# bio9.df <- raster::extract(climate$bio9,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio9, main = "Mean Temperature of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 10 ####
# bio10.df <- raster::extract(climate$bio10,             # raster layer
#                            points_spdf,   # SPDF with centroids for buffer
#                            buffer = 1,     # buffer size, units depend on CRS
#                            fun=mean,         # what to value to extract
#                            df=TRUE)         # return a dataframe? 

plot(climate$bio10, main = "Mean Temperature of Warmest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 11 ####
# bio11.df <- raster::extract(climate$bio11,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio11, main = "Mean Temperature of Coldest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 12 ####
# bio12.df <- raster::extract(climate$bio12,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio12, main = "Annual Precipitation") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 13 ####
# bio13.df <- raster::extract(climate$bio13,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio13, main = "Precipitation of Wettest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 14 ####
# bio14.df <- raster::extract(climate$bio14,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio14, main = "Precipitation of Driest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 15 ####
# bio15.df <- raster::extract(climate$bio15,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio15, main = "Precipitation Seasonality (Coefficient of Variation)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 16 ####
# bio16.df <- raster::extract(climate$bio16,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio16, main = "Precipitation of Wettest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 17 ####
# bio17.df <- raster::extract(climate$bio17,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio17, main = "Precipitation of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 18 ####
# bio18.df <- raster::extract(climate$bio18,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio18, main = "Precipitation of Warmest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 19 ####
# bio19.df <- raster::extract(climate$bio19,             # raster layer
#                             points_spdf,   # SPDF with centroids for buffer
#                             buffer = 1,     # buffer size, units depend on CRS
#                             fun=mean,         # what to value to extract
#                             df=TRUE)         # return a dataframe? 

plot(climate$bio19, main = "Precipitation of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)


climate.data <- cbind(bio1.df, bio2.df, bio3.df, bio4.df, bio5.df, bio6.df, bio7.df, bio8.df, bio9.df, bio10.df, bio11.df, bio12.df, bio13.df, bio14.df, bio15.df, bio16.df, bio17.df, bio18.df, bio19.df)

climate.data$culture_ID <- latlong$culture_ID
climate.data <- climate.data[-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)]

write.csv(climate.data, "output/climate_data.csv")
