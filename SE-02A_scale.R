## Scale code for lecture/lab
## SP25 Spatial Analysis Course
## 2024-08-01
## Erin Buchholtz, ekbuchh@clemson.edu
##------------------------------------------------------------------------------

## Notes & description ---------------------------------------------------------

# This code comes from Fletcher & Fortin textbook, "Spatial Ecology & 
# Conservation Modeling" second edition and updated for this course

# This code accompanies the lecture and serves to illustrate basic
# ideas of scale with raster data

## Libraries -------------------------------------------------------------------
library(terra)                  
library(sf)          
library(FedData)
library(tmap)
library(tidyverse)              

## -----------------------------------------------------------------------------

## Basics ----------------------------------------------------------------------
aea_proj <- "epsg:5070"
## -----------------------------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#2.3.3 A simple example
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(16)
toy <- rast(ncol=6, nrow=6, xmin=1, xmax=6, ymin=1, ymax=6)
toy[] <- rpois(ncell(toy), lambda=3)

#plot
plot(toy, axes=F)
text(toy, digits=2)

#check cell labeling/order
ncell(toy)
toy2 <- toy
toy2[] <- 1:ncell(toy)

#plot
plot(toy2, axes=F)
text(toy2, digits=2)

#increase the grain
toy_mean <- aggregate(toy, fact=2, fun="mean") #mean value
toy_maj <- aggregate(toy, fact=2, fun="modal") #majority rule

#plot mean rule
plot(toy_mean, axes=F)
text(toy_mean, digits=1)

#plot majority rule
plot(toy_maj, axes=F)
text(toy_maj)

#contrast means/variances
global(toy, mean) #or: mean(toy[])
global(toy, var)  #or: var(toy[])

global(toy_mean, mean)
global(toy_mean, var)

global(toy_maj, mean)
global(toy_maj, var)

#decrease the grain
toy_dis2 <- disagg(toy, fact=2)
toy_dis2_bilinear <- disagg(toy, fact=2, method='bilinear')

#plot
plot(toy_dis2, axes=F)
plot(as.polygons(toy_dis2, dissolve=F), add=TRUE, border='gray50', lwd=1)
text(toy_dis2, cex=0.9)

#plot
plot(toy_dis2_bilinear, axes=F)
plot(as.polygons(toy_dis2_bilinear, dissolve=F), add=TRUE, border='gray50', lwd=1)
text(toy_dis2_bilinear, digits=1, cex=0.6)

#decrease the extent
e <- ext(2, 4, 2, 4)     #first create new, smaller extent
toy_crop <- crop(toy, e) #then crop based on new extent

#plot
plot(toy)
plot(toy_crop)

#increase the extent
e <- ext(0, 7, 0, 7)      #first create new, bigger extent
toy_big <- extend(toy, e) #then crop based on new extent

#plot
plot(toy)
plot(toy_big)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#2.3.4.1 Multi-scale analysis - Prep
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#-------------------------------#
#site locations: shp file
#-------------------------------#

#site and reptile data
sites <- st_read("data/reptile_sites.shp") 

#inspect
class(sites)
crs(sites)
st_crs(sites) <- "epsg:5070" #set projection
crs(sites, describe = F, proj = T)
summary(sites)
head(sites, 2)

#study area including 10 km buffer around sampling points
studyarea <- st_as_sfc(st_bbox(sites) + c(-10000, -10000, 
                                            10000, 10000))
#view
tmap_mode("view")
tm_shape(studyarea)+tm_borders()+
  tm_shape(sites)+tm_dots()

#------------------#
#nlcd
#------------------#

# Download NLCD data for the defined area of interest
nlcd_p <- rast("data/nlcd2021_se.tif") |>
  project(aea_proj) 

writeRaster(nlcd, "data/nlcd2021_se_aea.tif")

nlcd <- nlcd %>%
  crop(studyarea) %>% mask(studyarea)

# nlcd <- get_nlcd(template = studyarea, 
#                       label = "NLCD", 
#                       year = 2021, 
#                       dataset = "landcover")
# 

# Print information about the downloaded raster
print(nlcd)

#inspect
crs(nlcd, proj=T)

#check same projection as reptile site data
crs(nlcd)==crs(sites)

#inspect raster properties
res(nlcd)
ncell(nlcd)
ext(nlcd)
levels(nlcd)

#reclassify to just land cover of interest



#reclassify with reclassify function is faster
levels(nlcd)
reclass_vals <- c(rep(0,7),  #NLCD classes 11-24
                  rep(1,3),  #NLCD classes 41-43
                  rep(0,10)) #NLCD classes 51-95
reclass_vals

#create reclassify matrix: first col=orginal value; second col=new class value
reclass_mat <- cbind(levels(nlcd)[[1]][,1], reclass)
reclass_mat

#reclassify nlcd to just 1/0 forest/not forest
forest <- classify(nlcd, reclass_mat)

plot(forest) #what is the current grain & extent?

#aggregate forest landcover to grain reflected in sampling
forest210 <- aggregate(forest, fact = 7, fun = modal)

#inspect
res(forest210) #210x210 res

#Write raster to file for next lab
writeRaster(forest210, "data/nlcd2016_se_forest210.tif")

#plot
plot(forest, axes = F)
plot(forest210, axes = F)
plot(sites, pch=20, col="blue", add=T)


#examine basic differences in buffers around sites
sitesamp_3km <- sites |> 
  filter(startsWith(site, "AL")) |>
  slice_sample(n = 10) |>
  st_buffer(dist = 3000) 

sitesamp_IDs <- sitesamp_3km[,1] |> st_drop_geometry()

# For the 30m res raster, view example sample site and calculate number of 
# pixels per land cover class within the site buffer polygons

s1_30 <- terra::crop(forest, sitesamp_3km[1,]) |>
  terra::mask(sitesamp_3km[1,])

sitesamp_forest30 <- terra::extract(forest, 
                                    sitesamp_3km,
                                    fun = table) |>
  rename(nonforest = "0", forest = "1") |>
  mutate(grain_m = 30) |>
  cbind(sitesamp_IDs)

# For the 200m res raster, view example sample site and calculate number of 
# pixels per land cover class within the site buffer polygons

s1_210 <- terra::crop(forest200, sitesamp_3km[1,]) |>
  terra::mask(sitesamp_3km[1,])

sitesamp_forest210 <- terra::extract(forest210, 
                                     sitesamp_3km,
                                     fun = table) |>
  rename(nonforest = "0", forest = "1") |>
  mutate(grain_m = 210) |>
  cbind(sitesamp_IDs)

# Join and pivot into long form, calculate area in km2
sitesamp_forestpix <- rbind(sitesamp_forest30, sitesamp_forest210) |>
  mutate(forest_perc = 100*round(forest/(nonforest+forest), digits = 3),
         tot_area_ha = ((nonforest+forest)*(grain_m*grain_m)/10000)) |>
  pivot_longer(cols = nonforest:forest, names_to = "cover", values_to = "count") |>
  mutate(area_ha = (grain_m*grain_m*count)/10000)

#sample differences
sitesamp_forestpix |> filter(ID == "1")
par(mfrow=c(1,2))
plot(s1_30)
plot(s1_210)
par(mfrow=c(1,1))

#plot differences
ggplot(sitesamp_forestpix, aes(x = site, y = forest_perc, color = as.factor(grain_m)))+
  geom_point(size = 4)+
  scale_color_manual(values = c("#578279", "#968C8C"))

ggplot(sitesamp_forestpix, aes(x = cover, y = area_ha))+
  geom_boxplot(aes(fill = as.factor(grain_m)), alpha = 0.5)+
  geom_jitter(aes(color = as.factor(grain_m)), position = position_dodge(width = 0.75), 
              alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("#578279", "#968C8C"))+
  scale_color_manual(values = c("black","black"))
