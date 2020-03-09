#reference https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html

pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, ggplot2, cowplot, tidyverse, rgdal, rasterVis)

#read the image file
base <- "../data/D3T13exp_Dec2019_Poorman/transect copy.tif"
cells <- "../data/D3T13exp_Dec2019_Poorman/cells.tif"
Fe <- "../data/D3T13exp_Dec2019_Poorman/Fe_transect.tif"
S <- "../data/D3T13exp_Dec2019_Poorman/S_transect.tif"

base_image <- image_read(base)

#convt to grayscale and normalize the intensity range 
Fe_gray <- image_read(Fe) %>% image_quantize(colorspace = 'gray') %>% image_equalize()
S_gray<- image_read(S) %>% image_quantize(colorspace = 'gray') %>% image_equalize()
cells_gray <- image_read(cells) %>% image_quantize(colorspace = 'gray') %>% image_equalize()


image_to_polygon <- function(image_path, to_raster = TRUE, drop_and_fill = TRUE, polygonize = TRUE){
  gray_image<- image_read(image_path) %>% image_quantize(colorspace = 'gray') %>% image_equalize()
  if(to_raster == TRUE && polygonize == FALSE && drop_and_fill == FALSE){
    temp_file <- tempfile()
    image_write(gray_image, path = temp_file, format = 'tiff')
    image_raster <- raster(temp_file) %>%
    cut(breaks = c(-Inf, 150, Inf)) - 1
    image_raster
  }else if(to_raster == TRUE && polygonize == TRUE && drop_and_fill == FALSE){
    image_polygon <- rasterToPolygons(image_raster, function(x){x == 1}, dissolve = TRUE) %>% 
      st_as_sf() 
    image_polygon
  }else if(to_raster == TRUE && polygonize == TRUE && drop_and_fill == TRUE){
    image_polygon %>% st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
    drop_crumbs(threshold = 100) %>%
    fill_holes(threshold = 100)
  }else{
    gray_image
  }
}

base_image_gray <- image_to_polygon(base, to_raster = F,drop_and_fill = F, polygonize = F)
Fe_raster <- image_to_polygon(Fe, drop_and_fill = F, polygonize = F)
S_raster <- image_to_polygon(S, drop_and_fill = F, polygonize = F)
cells_raster <- image_to_polygon(cells, drop_and_fill = F, polygonize = F)

plot(base_image_gray)
plot(Fe_raster,border = NA, lwd = 1.5, add = TRUE)

# polygonize
Fe_polygon <- rasterToPolygons(Fe_raster, function(x){x == 1}, dissolve = TRUE) %>% 
  st_as_sf() %>%
  st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
  drop_crumbs(threshold = 100) %>%
  fill_holes(threshold = 100)

# Fe_polygon <- st_as_sf(read_stars(Fe_file), as_points = TRUE, merge = FALSE)
Fe_polygon <- as(Fe_raster,'SpatialPolygonsDataFrame')
Fe_polygon <- fasterVectorize(Fe_raster, 'area', grassDir=grassDir)

S_polygon <- rasterToPolygons(S_raster, function(x){x == 1}, dissolve = TRUE) %>% 
  st_as_sf() %>%
  st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
  drop_crumbs(threshold = 100) %>%
  fill_holes(threshold = 100)


base_image_gray <- image_read(base) %>% image_quantize(colorspace = 'gray') %>% image_equalize()
base_temp_file <- tempfile()
image_write(base_image_gray, path = base_temp_file, format = 'tiff')
base_raster <- raster(base_temp_file)

plot(base_image_gray)
# plot(polygon, col= rgb(1, 0, 0,0.1), border = NA, lwd = 1.5, add = TRUE)
# plot(polygon2, col= rgb(0, 1, 1,0.1), border = NA, lwd = 1.5, add = TRUE)

pyrite <- st_intersection(st_union(polygon),st_union(polygon2)) #region of overlap between Fe and S 
S <- st_difference(st_union(polygon2),st_union(polygon)) #zone of just S
Fe <- st_difference(st_union(polygon),st_union(polygon2)) #zone of just Fe

plot(pyrite, col= rgb(1, 0, 1,0.1), border = NA, lwd = 1.5, add = TRUE)
plot(Fe, col= rgb(1, 0, 0,0.1), border = NA, lwd = 1.5, add = TRUE)
plot(S, col= rgb(1, 1, 0,0.1), border = NA, lwd = 1.5, add = TRUE)



# image stats -------------------------------------------------------------

cells_polygon <- rasterToPolygons(cells_raster, function(x){x == 1}, dissolve = TRUE) %>%
  st_as_sf() %>%
  st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
  fill_holes(threshold = 500) %>%
  drop_crumbs(threshold = 25) %>% #raster to polygon was creating some tiny polygons 
  st_cast("POLYGON")


# library(cluster)
# ellipsoid <- predict(ellipsoidhull(cells_polygon_coords))
# plot(cells_polygon_coords, asp=1)
# lines(ellipsoid)
# me <- colMeans((ellipsoid))   
# dist2center <- sqrt(rowSums((t(t(ellipsoid)-me))^2))
# max(dist2center)  
# min(dist2center)
# plot(ellipsoid,type='l',asp=1)
# points(cells_polygon_coords,col='blue')
# points(me,col='red')
# lines(rbind(me,ellipsoid[dist2center == min(dist2center),]))
# lines(ellipsoid[dist2center == max(dist2center),])

# 
# length <- tibble()
# width <- tibble()
# for (i in 1:(nrow(cells_polygon))) {
#   box <- unlist(unname(st_bbox(cells_polygon[i,])))
#   length <- append(length, box[3] - box[1])
#   width <- append(width, box[4] - box[2])
# }

#convert cell polygons to ellipses and calculate long + short axes
#this is imperfect due to curvature of longer cells 
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


ggplot(data = cells_polygon_stats) + 
  #geom_sf(aes(color = area)) +
  #scale_color_gradient(low = 'red', high = 'green') +
  geom_sf(aes(color = shape))

# calculate the centroid 
cell_centroids <- st_centroid(cells_polygon) %>%
  mutate(geom = str_replace_all(geometry, "c|[)]|[(]", "")) %>%
  separate(geom, c("x", "y"), sep = ",") %>%
  mutate(x = as.numeric(x), y = as.numeric(y),
         shape = cells_polygon_stats$shape,
         area = cells_polygon_stats$area) %>%
  as.tibble() %>%
  select(x, y, shape, area) #%>%
  # mutate(z = 1) %>%
  # as.matrix() %>%
  # SpatialPoints()
  

plot(base_image)
points(y~x, data = sf_cent, pch=20,col = "red")

##plot with gplot 

base_image <- raster(base) 
Fe_raster_coords <- as.data.frame(Fe_raster, xy=TRUE)
S_raster_coords <- as.data.frame(S_raster, xy=TRUE)
rasterVis::gplot(base_image) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  ggnewscale::new_scale_fill() +
  geom_raster(Fe_raster_coords, mapping = aes(x, y, fill = layer), alpha = 0.3) +
  scale_fill_gradient(low = NA, high = 'red') +
  ggnewscale::new_scale_fill() +
  geom_raster(S_raster_coords, mapping = aes(x, y, fill = layer), alpha = 0.3) +
  scale_fill_gradient(low = NA, high = 'yellow') +
  coord_fixed() +
  #geom_point(data=cell_centroids, mapping = aes(x, y), size = .1, col="orange") +
  geom_sf(data = cells_polygon, inherit.aes = F, color = "orange") +
  ggnewscale::new_scale_color() +
  geom_density_2d(data=cell_centroids, inherit.aes = F, mapping = aes(x, y, col = stat(level)/max(stat(level)))) +
  scale_color_viridis_c() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())


##plot only with ggplot 
base_coords <- as.data.frame(base_image, xy=TRUE)
Fe_raster_coords <- as.data.frame(Fe_raster, xy=TRUE)
ggplot() +
  geom_tile(Fe_raster_coords, mapping = aes(x, y, fill = layer), color = "red") +
  coord_fixed()

###plot with raserviz

# Set color palette
zeroCol <-NA 

FeTheme <- rasterTheme(region = c(zeroCol, 'red'))
STheme <- rasterTheme(region = c(zeroCol, 'yellow'))
cellTheme <- rasterTheme(region = c(zeroCol, 'blue'))



background_image <- levelplot(base_image, par.settings = GrTheme, margin = FALSE)
Fe_map <- levelplot(Fe_raster, par.settings = FeTheme, margin = FALSE, alpha.regions = 0.35)
S_map <- levelplot(S_raster, par.settings = STheme, margin = FALSE, alpha.regions = 0.35)
cell_map <- levelplot(cells_raster, par.settings = cellTheme, margin = FALSE, alpha.regions = 0.35, contour = TRUE)

background_image + Fe_map + S_map + cell_map +
  layer(sp.points(cell_centroids, pch=20, cex=0.1, col=1), columns=1) 
  

densityplot(cells_raster)



# stitch images -----------------------------------------------------------



library(image.OpenPano)

folder <- system.file(package = "image.OpenPano", "extdata")
images <- c(file.path(folder, "imga.jpg"), 
            file.path(folder, "imgb.jpg"),
            file.path(folder, "imgc.jpg"))
result <- image_stitch(images, file = "result_stitched.jpg")

library(magick)
image_read(images[1])
image_read(images[2])
image_read(images[3])
image_read("result_stitched.jpg")




gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                  pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}




# spatial stats -----------------------------------------------------------
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

#cells centroids
cell_centroids <- st_centroid(cells_polygon)
cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame()

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


