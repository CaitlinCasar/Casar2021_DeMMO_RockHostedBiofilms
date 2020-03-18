pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, ggplot2, cowplot, tidyverse, rgdal, rasterVis)

experiment_id <- "DeMMO3/D3T13exp_Dec2019_Poorman"

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

#store images in appropriate format for analyses
base_image <- image_to_polygon(base, polygonize = F, pres_abs = F)
Fe_raster <- image_to_polygon(Fe, polygonize = F) #should return intensities between 0-1
S_raster <- image_to_polygon(S, polygonize = F)
cells_polygon <- image_to_polygon(cells) 

S_polygon <- rasterToPolygons(S_raster, function(x){x == 1}, dissolve = TRUE) %>% 
  st_as_sf() %>%
  st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
  drop_crumbs(threshold = 100) %>%
  fill_holes(threshold = 100)

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
  #ggnewscale::new_scale_fill() +
  #geom_raster(as.data.frame(Fe_raster, xy=TRUE), mapping = aes(x, y, fill = layer), alpha = 0.3) +
  #scale_fill_gradient(low = NA, high = 'red') +
  ggnewscale::new_scale_fill() +
  geom_raster(as.data.frame(S_raster, xy=TRUE), mapping = aes(x, y, fill = layer), alpha = 0.3) +
  scale_fill_gradient(low = NA, high = 'yellow') +
  coord_fixed() +
  geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
  ggnewscale::new_scale_color() +
  geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
  scale_color_viridis_c() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())

# spatial stats -----------------------------------------------------------
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

#create ppp for point pattern analysis in spatstat
window <- owin(xrange = c(st_bbox(cell_centroids)$xmin, st_bbox(cell_centroids)$xmax),
               yrange = c(st_bbox(cell_centroids)$ymin, st_bbox(cell_centroids)$ymax))
cell_centroids_ppp <- ppp(cell_centroids_coords$X, cell_centroids_coords$Y, window) %>%
  rescale(22.15, "μm")

#compute quadrat density
cells_quadrat <- quadratcount(cell_centroids_ppp, nx= 12, ny=3)
plot(cell_centroids_ppp, main=NULL, cols=rgb(0,0,0,.2), pch=20)
plot(cell_centroids_ppp, pch=20, cols="grey70", main=NULL)  # Plot points
plot(cells_quadrat, add=TRUE)  # Add quadrat grid
# Compute the density for each quadrat
cells_quadrat_dens <- intensity(cells_quadrat)

# Plot the density
plot(intensity(cells_quadrat, image=TRUE), main=NULL, las=1)  # Plot density raster
plot(cell_centroids_ppp, pch=20, cex=0.1, col=rgb(0,0,0,.5), add=TRUE)  # Add points


cell_density <- density(cell_centroids_ppp) # Using the default bandwidth
plot(cell_density, main=NULL, las=1)
contour(cell_density, add=TRUE)


##point process model
#test whether elements explain cell densities
S_im <- as.im(S_raster) %>%
  rescale(22.15, "μm")

Fe_im <- as.im(Fe_raster) %>%
  rescale(22.15, "μm")

point_process_model1 <- ppm(cell_centroids_ppp ~ S_im + Fe_im)
point_process_model0 <- ppm(cell_centroids_ppp ~ 1)
anova(point_process_model0, point_process_model1, test="LRT")


#calculate spatial randomness 
#do we expect to see clustering?
#Where K falls under the theoretical Kpois line the points are more clustered at distance r, and vis versa

cell_kest <- Kest(cell_centroids_ppp)
plot(cell_kest, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))

ANN <- apply(nndist(cell_centroids_ppp, k=1:100),2,FUN=mean)
plot(ANN ~ eval(1:100), type="b", main=NULL, las=1)

ann.p <- mean(nndist(cell_centroids_ppp, k=1))
ann.p


#generate random distribution of cells as null model
n     <- 599L               # Number of simulations
ann.r <- vector(length = n) # Create an empty object to be used to store simulated ANN values
for (i in 1:n){
  rand.p   <- rpoint(n=cell_centroids_ppp$n, win=window)  # Generate random point locations
  ann.r[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}

plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))

hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

#test whether S distribution explains cell distribution 
n     <- 599L
ann.r <- vector(length=n)
for (i in 1:n){
  rand.p   <- rpoint(n=cell_centroids_ppp$n, f=S_im) 
  ann.r[i] <- mean(nndist(rand.p, k=1))
}

Window(rand.p) <- window  # Replace raster mask with ma.km window
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))

hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p





# calculate % element per pixel -------------------------------------------
#need to automate x-ray image stitching for each transect
Fe_intensity <- raster(Fe)
S_intensity <- raster(S)

SF_raster <- stack(c(Fe_intensity, S_intensity))

Fe_percent <- overlay(SF_raster, fun=function(Fe, S)Fe/(Fe+S))
S_percent <- overlay(SF_raster, fun=function(Fe, S)S/(Fe+S))

element_perc_df <- as.data.frame(rasterToPoints(Fe_intensity, xy = TRUE)) %>%
  left_join(as.data.frame(rasterToPoints(S_intensity, xy = TRUE)))

test <- element_perc_df %>%
  mutate(Fe_S = Fe/S) %>%
  #filter(Fe_S > 0.45 & Fe_S < 0.55) %>%
  dplyr::select(-Fe, -S) %>%
  rasterFromXYZ() 

test1 <- test %>%
  cut(breaks = c(0.45, 0.5, 0.55)) -1

plot(test)

pyrite <- rasterToPolygons(test1, function(x){x == 1}, dissolve = T) %>% 
  st_as_sf() 


if(drop_and_fill == TRUE){
  image_polygon_drop_fill <- image_polygon %>% 
    st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
    fill_holes(threshold = 500) %>%
    drop_crumbs(threshold = 25) %>% #raster to polygon was creating some tiny polygons 
    st_cast("POLYGON")
  
  