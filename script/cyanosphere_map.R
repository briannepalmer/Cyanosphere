# make map of cyanosphere points

# load libraries 
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(tidyverse)
library(ggrepel)
library(paletteer)
library(viridis)

cyanosphere <- read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/Cyanosphere/data/metdata_NEW.csv")
fitlered.cyano <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
fitlered.cyano <- fitlered.cyano %>% filter(Phylum != "p__Cyanobacteria", Phylum != "p__Cyanobacteriota")

cyanosphere <- cyanosphere %>% filter(Host %in% fitlered.cyano$Host)
cyanosphere$Label <- c(1:56)

pal <- paletteer_d("ggsci::nrc_npg")
pal1 <- colorRampPalette(pal)(20)

# Group by culture_ID to get just 50 points 

world <- map_data("world") 
gg1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "gray") + 
  coord_fixed(1.3) + theme_bw()

png("figures/map.png", res = 300, units = "cm", height = 15, width = 24)
gg1 + 
  geom_point(data = cyanosphere, aes(x = Long, y = Lat, color = Substrate), size = 2) + geom_jitter() + geom_text_repel(data = cyanosphere, aes(x = Long, y = Lat, label = Label, color = Substrate), max.overlaps = 100) + scale_color_manual(values = pal) + theme(text=element_text(size=20))
dev.off()


# add climate data
library(raster)
library(geodata)

# run to get climate data
climate <- worldclim_global(var = 'bio', res = 2.5, path=tempdir())

#### Bio 1 ####
latlong <- cyanosphere %>% group_by(Long, Lat, Host) %>% summarize()

latlong <- na.omit(latlong)

point <- SpatialPoints(coords = cbind(latlong$Long,latlong$Lat))
crs <- crs(point)

points_spdf <- SpatialPointsDataFrame(latlong[,1:2], proj4string= crs, latlong)

bio1.df <- raster::extract(climate$bio1,points_spdf,buffer = 1,fun=mean,  df=TRUE)         # return a dataframe?

# mean temp is C degrees times ten

bio1.df$bio1 <- bio1.df$bio1/10

plot(climate$wc2.1_2.5m_bio_1, main = "Annual Mean Temperature") 
plot(climate$wc2.1_2.5m_bio_1, main = "Annual Mean Temperature", xlim = c(- 125, -100), ylim = c(30, 40)) 
plot(climate$wc2.1_2.5m_bio_12, main = "Annual Precipitation",, xlim = c(- 125, -100), ylim = c(30, 40)) 

test <- raster(climate$wc2.1_2.5m_bio_1)
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df1 <- as.data.frame(test_spdf)
colnames(test_df1) <- c("value", "x", "y")

us <- map_data("state")
us.west <- us %>% filter(region == c("california", "nevada", "arizona", "new mexico", "utah", "colorado"))
test_df1 <- test_df1 %>% filter(y < max(us.west$lat), x < max(us.west$long), y > min(us.west$lat), x > min(us.west$long))


ggplot() +
  geom_polygon(
    data = us,
    aes(x = long, y = lat, group = group),
    fill = "white", color = "black"
  )  + scale_fill_viridis() + theme_bw() + coord_map("albers", lat0=32.5343, lat=142.0095) + geom_tile(data = test_df1, aes(x = x, y = y, fill = value), alpha = 0.8)

test <- raster(climate$wc2.1_2.5m_bio_12)
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df2 <- as.data.frame(test_spdf)
colnames(test_df2) <- c("value", "x", "y")
test_df2 <- test_df2 %>% filter(y < max(us.west$lat), x < max(us.west$long), y > min(us.west$lat), x > min(us.west$long))

ggplot() +
  geom_polygon(
    data = us,
    aes(x = long, y = lat, group = group),
    fill = "white", color = "black"
  )  + scale_fill_viridis(option = "A") + theme_bw() + coord_map("albers", lat0=32.5343, lat=142.0095) + geom_tile(data = test_df2, aes(x = x, y = y, fill = value), alpha = 0.8)

ggplot() +
  geom_polygon(
    data = us,
    aes(x = long, y = lat, group = group),
    fill = "white", color = "black"
  )  + scale_fill_viridis(option = "A") + theme_bw() + coord_map("albers", lat0=32.5343, lat=142.0095) + geom_tile(data = test_df1, aes(x = x, y = y, fill = value), alpha = 0.8) 



points(x = latlong$Long, y = latlong$Lat)


#### Bio 2 ####
bio2.df <- raster::extract(climate$bio2,             # raster layer
                               points_spdf,   # SPDF with centroids for buffer
                               buffer = 1,     # buffer size, units depend on CRS
                               fun=mean,         # what to value to extract
                               df=TRUE)         # return a dataframe?


plot(climate$bio2, main = "Mean Diurnal Range (Mean of monthly (max temp - min temp))") 


points(x = latlong$Long, y = latlong$Lat)


#### Bio 3 ####
bio3.df <- raster::extract(climate$bio3,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio3, main = "Isothermality (BIO2/BIO7) (×100)")
points(x = latlong$Long, y = latlong$Lat)

#### Bio 4 ####
bio4.df <- raster::extract(climate$bio4,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio4, main = "Temperature Seasonality (standard deviation ×100)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 5 ####
bio5.df <- raster::extract(climate$bio5,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio5, main = "Max Temperature of Warmest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 6 ####
bio6.df <- raster::extract(climate$bio6,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio6, main = "Min Temperature of Coldest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 7 ####
bio7.df <- raster::extract(climate$bio7,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio7, main = "Temperature Annual Range (BIO5-BIO6)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 8 ####
bio8.df <- raster::extract(climate$bio8,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio8, main = "Mean Temperature of Wettest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 9 ####
bio9.df <- raster::extract(climate$bio9,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio9, main = "Mean Temperature of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 10 ####
bio10.df <- raster::extract(climate$bio10,             # raster layer
                           points_spdf,   # SPDF with centroids for buffer
                           buffer = 1,     # buffer size, units depend on CRS
                           fun=mean,         # what to value to extract
                           df=TRUE)         # return a dataframe?

plot(climate$bio10, main = "Mean Temperature of Warmest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 11 ####
bio11.df <- raster::extract(climate$bio11,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio11, main = "Mean Temperature of Coldest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 12 ####
bio12.df <- raster::extract(climate$bio12,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio12, main = "Annual Precipitation") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 13 ####
bio13.df <- raster::extract(climate$bio13,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio13, main = "Precipitation of Wettest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 14 ####
bio14.df <- raster::extract(climate$bio14,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio14, main = "Precipitation of Driest Month") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 15 ####
bio15.df <- raster::extract(climate$bio15,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio15, main = "Precipitation Seasonality (Coefficient of Variation)") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 16 ####
bio16.df <- raster::extract(climate$bio16,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio16, main = "Precipitation of Wettest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 17 ####
bio17.df <- raster::extract(climate$bio17,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio17, main = "Precipitation of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 18 ####
bio18.df <- raster::extract(climate$bio18,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio18, main = "Precipitation of Warmest Quarter") 
points(x = latlong$Long, y = latlong$Lat)

#### Bio 19 ####
bio19.df <- raster::extract(climate$bio19,             # raster layer
                            points_spdf,   # SPDF with centroids for buffer
                            buffer = 1,     # buffer size, units depend on CRS
                            fun=mean,         # what to value to extract
                            df=TRUE)         # return a dataframe?

plot(climate$bio19, main = "Precipitation of Driest Quarter") 
points(x = latlong$Long, y = latlong$Lat)


climate.data <- cbind(bio1.df, bio2.df, bio3.df, bio4.df, bio5.df, bio6.df, bio7.df, bio8.df, bio9.df, bio10.df, bio11.df, bio12.df, bio13.df, bio14.df, bio15.df, bio16.df, bio17.df, bio18.df, bio19.df)

climate.data$Host <- latlong$Host
climate.data <- climate.data[-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)]

write.csv(climate.data, "output/climate_data.csv")

# merge with metdatdata 

metadata <- rbind(cyanosphere, climate.data,)
metadata <- merge(cyanosphere, climate.data, by = "Host", all = TRUE)

write.csv(metadata, "data/metdata_Feb2024.csv")



########################

cyanosphere <- read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/Cyanosphere/data/metdata_NEW.csv")
fitlered.cyano <- read.csv("data/cyanosphere_clean_2024_lowCheckM_removed.csv")
fitlered.cyano <- fitlered.cyano %>% filter(Phylum != "p__Cyanobacteria", Phylum != "p__Cyanobacteriota", Habitat2 == "desert_soil")

cyanosphere <- cyanosphere %>% filter(Host %in% fitlered.cyano$Host)
cyanosphere$Label <- c(1:24)

pal <- paletteer_d("ggsci::nrc_npg")
pal1 <- colorRampPalette(pal)(20)

# Group by culture_ID to get just 50 points 

world <- map_data("world") 
gg1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "gray") + 
  coord_fixed(1.3) + theme_bw()

gg1 + 
  geom_point(data = cyanosphere, aes(x = Long, y = Lat, color = Location), size = 2) + geom_jitter() + geom_text_repel(data = cyanosphere, aes(x = Long, y = Lat, label = Label, color = Location), max.overlaps = 100) + scale_color_manual(values = pal) + theme(text=element_text(size=20))



gg1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "gray") + 
  coord_fixed(1.3) + theme_bw()

gg1 + 
  geom_point(data = cyanosphere, aes(x = Long, y = Lat, color = Location), size = 2) + geom_jitter() + geom_text_repel(data = cyanosphere, aes(x = Long, y = Lat, label = Label, color = Location), max.overlaps = 100) + scale_color_manual(values = pal) + theme(text=element_text(size=20))

## make map with aridity ##

#--- Download the files from the TerraClimate website ---#
# Precipitation
download.file(url = 'http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_ppt_2019.nc',
              destfile = 'ppt.nc')

# Evapotranspiration
download.file(url = 'http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_pet_2019.nc',
              destfile = 'pet.nc')

#--- Import the downloaded files ---#
# Precipitation
ppt <- stack(x = 'ppt.nc')

# Evapotranspiration
pet <- stack(x = 'pet.nc')

#--- Inspect ---#
# Precipitation
plot(ppt)

#--- Raster maths ---#
# Precipitation
ppt_mean <- calc(ppt, # RasterStack object
                 fun = mean, # Function to apply across the layers
                 na.rm = TRUE)

# Evapotranspiration
pet_mean <- calc(pet,
                 fun = mean, 
                 na.rm = TRUE)

#--- Set the extent ---#
# Cut off all values below 60 degrees South (removing Antarctica)
ext <- extent(c(xmin = -180, xmax = 180, 
                ymin = -60, ymax = 90))

#--- Crop ---#
# Precipitation
ppt_mean <- crop(x = ppt_mean, 
                 y = ext)

# Evapotranspiration
pet_mean <- crop(x = pet_mean, 
                 y = ext)

#--- Inspect ---#
# Precipitation
plot(main = 'Precipitation',
     ppt_mean)

# Evapotranspiration
plot(main = 'Evapotranspiration',
     pet_mean)

#--- Calculate aridity index ---#
# Precipitation (ppt) / Evapotranspiration (pet)
aridity_index <- overlay(x = ppt_mean, # Raster object 1
                         y = pet_mean, # Raster object 2
                         fun = function(x, y){return(x / y)}) # Function to apply

#--- Convert raster to a matrix ---#
aridity_index_matrix <- rasterToPoints(aridity_index)

#--- Convert to the matrix to a dataframe ---#
aridity_index_df <- as.data.frame(aridity_index_matrix)

#--- Recode aridity index into categories --#
aridity_index_df <- aridity_index_df %>% 
  # Recode
  mutate(category = case_when(
    is.infinite(layer) ~ 'Humid',
    layer >= 0.65 ~ 'Humid',
    layer >= 0.5 & layer < 0.65 ~ 'Dry sub-humid',
    layer >= 0.2 & layer < 0.5 ~ 'Semi-arid',
    layer >= 0.05 & layer < 0.2 ~ 'Arid',
    layer < 0.05 ~ 'Hyper-arid'
  )) %>% 
  # Convert to ordered factor
  mutate(category = factor(category,
                           levels = c('Hyper-arid', 'Arid', 'Semi-arid',
                                      'Dry sub-humid', 'Humid'),
                           ordered = TRUE))

#--- Set a colour palette ---#
colours <- c('#e31a1c', '#fd8d3c', '#fecc5c', '#ffffb2', '#666666')

#--- Plot the data ---#
map <- ggplot(data = aridity_index_df) +
  aes(y = y,
      x = x,
      fill = category) +
  geom_raster() +
  scale_fill_manual(values = colours,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(limits = c(-60, 90),
                     expand = c(0, 0),
                     breaks = c(-40, -20, 0, 20, 40, 60, 80),
                     labels = c(expression('40'*degree*'S'),
                                expression('20'*degree*'S'),
                                expression('0'*degree),
                                expression('20'*degree*'N'),
                                expression('40'*degree*'N'),
                                expression('60'*degree*'N'),
                                expression('80'*degree*'N'))) +
  scale_x_continuous(limits = c(-180, 180),
                     expand = c(0, 0),
                     breaks = c(-180, -120, -60, 0, 60, 120, 180),
                     labels = c(expression('180'*degree*'W'),
                                expression('120'*degree*'W'),
                                expression('60'*degree*'W'),
                                expression('0'*degree),
                                expression('60'*degree*'E'),
                                expression('120'*degree*'E'),
                                expression('180'*degree*'E'))) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid.major = element_line(linetype = 2, 
                                        size = 0.5,
                                        colour = '#666666'),
        panel.grid.minor = element_blank()) 


ggsave("aridity_map.svg", map, device = "svg", height = 10, width = 15, units = "cm", dpi = 300)




