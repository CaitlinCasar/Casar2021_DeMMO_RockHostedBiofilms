pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, ggplot2, cowplot, tidyverse, rgdal, rasterVis)

experiment_id <- "D3T13exp_Dec2019_Poorman"

#assign image file paths 
base <- paste0("../data/", experiment_id, "/transect.tif")
cells <- paste0("../data/", experiment_id, "/cells.tif")
Fe <- paste0("../data/", experiment_id, "/Fe.tif")
S <- paste0("../data/", experiment_id, "/S.tif")


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
        image_polygon <- rasterToPolygons(image_raster, function(x){x == 1}, dissolve = TRUE) %>% 
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

#store images in appropriate format for analyses
base_image <- image_to_polygon(base, polygonize = F, pres_abs = F)
Fe_raster <- image_to_polygon(Fe, polygonize = F) #this function gives different output for Fe for some reason...  
S_raster <- image_to_polygon(S, polygonize = F)
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

#plot data with cell kernel density contours 
base_raster <- raster(base) 

density_plot <- rasterVis::gplot(base_raster) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  ggnewscale::new_scale_fill() +
  geom_raster(as.data.frame(Fe_raster, xy=TRUE), mapping = aes(x, y, fill = layer), alpha = 0.3) +
  scale_fill_gradient(low = NA, high = 'red') +
  ggnewscale::new_scale_fill() +
  geom_raster(as.data.frame(S_raster, xy=TRUE), mapping = aes(x, y, fill = layer), alpha = 0.3) +
  scale_fill_gradient(low = NA, high = 'yellow') +
  coord_fixed() +
  geom_sf(data = cells_polygon, inherit.aes = F, color = "orange") +
  ggnewscale::new_scale_color() +
  geom_density_2d(data=cell_centroids, inherit.aes = F, mapping = aes(x, y, col = stat(level)/max(stat(level)))) +
  scale_color_viridis_c() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())

# spatial stats -----------------------------------------------------------
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html




