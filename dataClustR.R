#load dependencies
pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, ggplot2, cowplot, tidyverse, rgdal, rasterVis)

#set working directory
setwd("/Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO1/D1T1exp_Dec2019_Poorman")

#load xray raster brick 
xray_brick <- brick("~/Desktop/dataStitcher/example_dataset/example_dataset_brick.grd")

#cell image path
cells <- "cells.tif"

#function for rasterizing or polygonizing
image_to_polygon <- function(image_path, to_raster = TRUE, pres_abs = TRUE, drop_and_fill = TRUE, polygonize = TRUE, equalizer = TRUE){
  if(equalizer == T){
    image <- image_path %>% image_read() %>% image_quantize(colorspace = 'gray') %>% image_equalize() 
    image
    if(to_raster == T){
      temp_file <- tempfile()
      image_write(image, path = temp_file, format = 'tiff')
      image_raster <- raster(temp_file) %>%
        cut(breaks = c(-Inf, 150, Inf)) - 1
      image_raster 
      if(polygonize == TRUE){
        image_polygon <- rasterToPolygons(image_raster, function(x){x == 0}, dissolve = TRUE) %>% 
          st_as_sf() 
        if(drop_and_fill == TRUE){
          image_polygon_drop_fill <- image_polygon %>% 
            st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
            fill_holes(threshold = 500) %>%
            drop_crumbs(threshold = 25) %>% #raster to polygon was creating some tiny polygons 
            st_cast("POLYGON")
          image_polygon_drop_fill
        }else{
          image_polygon
        }
      }else{
        image_raster
      }
    }else{
      image
    }
  }else{
    image_read(image_path)
  }
}

#generate cell polygons
cells_polygon <- image_to_polygon(cells) 

#function to draw ellipses around cell polygons, measure major and minor axes 
ellipsoid_fun <- function(data){
  polygon_coords <- data %>%
    st_coordinates() %>%
    as.data.frame() %>%
    select(-L1, -L2) %>%
    as.matrix() 
  ellipsoid <- predict(ellipsoidhull(polygon_coords))
  me <- colMeans((ellipsoid))   
  dist2center <- sqrt(rowSums((t(t(ellipsoid)-me))^2))
  long_axis <- max(dist2center)  
  short_axis <- min(dist2center)
  long_axis/short_axis
}

#generate table of cell stats
cells_polygon_stats <- cells_polygon %>%
  mutate(area = st_area(geometry)/(22.15^2),
         perimeter = st_perimeter(geometry)/22.15)
cells_polygon_stats$roundness <- sapply(cells_polygon$geometry, ellipsoid_fun)
cells_polygon_stats <- cells_polygon_stats %>%
  mutate(shape = case_when(roundness < 1.12 ~ "cocci",
                           roundness >= 1.12 & roundness < 5 ~"rod",
                           TRUE ~ "filament"))

cells_polygon_totals <- cells_polygon_stats %>% 
  group_by(shape) %>% 
  summarise(n())

# calculate cell centroids 
cell_centroids <- st_centroid(cells_polygon)

cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(shape = cells_polygon_stats$shape,
         area = cells_polygon_stats$area) 


      