#load dependencies
pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, ggplot2, cowplot, tidyverse, rgdal, rasterVis)

#set working directory
setwd("/Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO1/D1T1exp_Dec2019_Poorman")

#load xray raster brick 
#xray_brick_not_equalized_or_thresholded <- brick("not_equalized/example_brick.grd") #wrong, Si showing zero
#xray_brick_equalized_not_thresholded <- brick("equalized_not_thresholded/example_brick.grd") #same as thresholded 
xray_brick <- brick("equalized_and_thresholded/example_brick.grd") #7% SI

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
            #st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
            st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +units=m") #need to check that this works, otherwise reproject 
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

cells_polygon  <- cells_polygon %>% st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# cell_image <- cells %>%
#   image_read() %>% 
#   image_quantize(colorspace = 'gray') %>% 
#   image_equalize() 
# temp_file <- tempfile()
# image_write(cell_image, path = temp_file, format = 'tiff')
# cell_raster <- raster(temp_file) %>%
#   cut(breaks = c(-Inf, 150, Inf)) - 1

#check plotting
# cells_filtered <- as.data.frame(cells_polygon, xy=TRUE) %>%
#   filter(layer == 0)
# rasterVis::gplot(flip(xray_brick[[1]], direction = 'y')) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = 'black', high = 'white') +
#   coord_fixed() +
#   geom_tile(cells_filtered, mapping = aes(x, y), fill = "gold")


rasterVis::gplot(xray_brick[[1]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  coord_fixed() +
  geom_sf(data = cells_polygon_projected, inherit.aes = F, fill = "orange", lwd = 0) +
  coord_sf(datum = NA)

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
  as.data.frame %>%
  group_by(shape) %>% 
  summarise(value = n()) %>%
  rename(observation = shape) %>%
  bind_rows(data.frame(observation = "mean cell area (μm)", value = mean(cells_polygon_stats$area)),
            data.frame(observation = "total cells", value = nrow(cells_polygon_stats)),
            data.frame(observation = "cell density (cells/mm^2)", value = nrow(cells_polygon_stats)/((240*33.66516)/1000000)))

# calculate cell centroids 
cell_centroids <- st_centroid(cells_polygon)

cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(shape = cells_polygon_stats$shape,
         area = cells_polygon_stats$area) 

# plot the data -----------------------------------------------------------
xy <- read_delim("coordinates.txt", delim = "\t")

SEM_images <- list.files("/Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO1/D1T1exp_Dec2019_Poorman/SEM_images/", full.names = T)
SEM_images <- SEM_images[!str_detect(SEM_images, "overview")]
SEM_panorama <- list()


message("Stitching SEM images into panorama...")
for(i in 1:length(SEM_images)){
  image <- SEM_images[i] %>% image_read() %>%
    image_quantize(colorspace = 'gray') %>%
    image_equalize()
  temp_file <- tempfile()
  image_write(image, path = temp_file, format = 'tiff')
  image <- raster(temp_file)
  image_extent <- extent(matrix(c(xy$x[i], xy$x[i] + 1024, xy$y[i], xy$y[i]+704), nrow = 2, ncol = 2, byrow = T))
  image_raster <- setExtent(raster(nrows = 704, ncols = 1024), image_extent, keepres = F)
  values(image_raster) <- values(image)
  SEM_panorama[[i]] <- image_raster
}
SEM_panorama_merged <- do.call(merge, SEM_panorama)
message("...complete.")




message("Generating plot...")
# Set color palette

#palette source: https://sciencenotes.org/molecule-atom-colors-cpk-colors/

element_colors <- c("#FFFFFF", "#D9FFFF", "#CC80FF", "#C2FF00", "#FFB5B5", "#909090", "#3050F8",
                    "#FF0D0D", "#90E050", "#B3E3F5", "#AB5CF2", "#8AFF00", "#BFA6A6", "#F0C8A0",
                    "#FF8000", "#FFFF30", "#1FF01F", "#80D1E3", "#8F40D4", "#3DFF00", "#E6E6E6",
                    "#BFC2C7", "#A6A6AB", "#8A99C7", "#9C7AC7", "#E06633", "#F090A0", "#50D050",
                    "#C88033", "#7D80B0", "#C28F8F", "#668F8F", "#BD80E3", "#FFA100", "#A62929",
                    "#5CB8D1", "#702EB0", "#00FF00", "#94FFFF", "#94E0E0", "#73C2C9", "#54B5B5",
                    "#3B9E9E", "#248F8F", "#0A7D8C", "#006985", "#C0C0C0", "#FFD98F", "#A67573",
                    "#668080", "#9E63B5", "#D47A00", "#940094", "#429EB0", "#57178F", "#00C900",
                    "#70D4FF", "#FFFFC7", "#D9FFC7", "#C7FFC7", "#A3FFC7", "#8FFFC7", "#61FFC7",
                    "#45FFC7", "#30FFC7", "#1FFFC7", "#00FF9C", "#00E675", "#00D452", "#00BF38",
                    "#00AB24", "#4DC2FF", "#4DA6FF", "#2194D6", "#267DAB", "#266696", "#175487",
                    "#D0D0E0", "#FFD123", "#B8B8D0", "#A6544D", "#575961", "#9E4FB5", "#AB5C00",
                    "#754F45", "#428296", "#420066", "#007D00", "#70ABFA", "#00BAFF", "#00A1FF",
                    "#008FFF", "#0080FF", "#006BFF", "#545CF2", "#785CE3", "#8A4FE3", "#A136D4",
                    "#B31FD4", "#B31FBA", "#B30DA6", "#BD0D87", "#C70066", "#CC0059", "#D1004F",
                    "#D90045", "#E00038", "#E6002E", "#EB0026")
names(element_colors) <- c("H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
                           "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
                           "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo",
                           "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
                           "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                           "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
                           "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                           "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt")

xray_frame <- as.data.frame(xray_brick, xy=TRUE) 
xray_frame <- gather(xray_frame, element, value, colnames(xray_frame)[3]:colnames(xray_frame)[ncol(xray_frame)])


element_plotter<-function(coord_frame, brick, SEM_image, colors){
  p <-rasterVis::gplot(SEM_image) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'black', high = 'white') +
    ggnewscale::new_scale_fill()
  for(i in names(brick)){ 
    print(paste0("Adding ", names(brick[[i]]), " to plot..."))
    element_coords <- coord_frame %>%
      filter(element == names(brick[[i]]) & value!=0)
    p <- p+geom_raster(element_coords, mapping = aes(x, y, fill = element, alpha = value)) +
      scale_fill_manual(values = colors) +
      ggnewscale::new_scale_fill()
  }
  message("Writing plot...")
  suppressWarnings(print(p + 
                           coord_fixed() +
                           geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
                           ggnewscale::new_scale_color() +
                           geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
                           scale_color_viridis_c() +
                           theme(axis.title = element_blank(),
                                 axis.text = element_blank(),
                                 legend.position = "none")))
}

element_plot <- element_plotter(xray_frame, xray_brick, SEM_panorama_merged, element_colors)
element_plot

element_plot_legend <- data.frame(element = unique(xray_frame$element)) %>%
  rownames_to_column() %>% 
  ggplot(aes(element, rowname, fill=element)) + 
  geom_bar(stat= "identity") + 
  scale_fill_manual(values = element_colors) +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8))

element_plot_with_legend <- plot_grid(
  element_plot, 
  plot_grid(get_legend(element_plot_legend), 
            ncol = 1), 
  nrow = 2, 
  rel_heights = c(8,2)
)


# generate PDF ------------------------------------------------------------
sample_id <- "D1T1exp_Dec2019_Poorman"
pdf(paste0(sample_id, "_density_plot.pdf"),
    width = 13.33, 
    height = 7.5)

print(element_plot_with_legend)

dev.off()

message("...complete.")

#remove(list = ls())


# summarize rock composition ----------------------------------------------
message("Generating reports...")

#per pixel summary 
# xray_brick_summary <- as.data.frame(xray_brick, xy=T) %>%
#   select(-x, -y)
# xray_brick_summary[is.na(xray_brick_summary)] <- 0
# 
# xray_brick_summary <- xray_brick_summary/rowSums(xray_brick_summary) * 100
# xray_brick_summary[is.na(xray_brick_summary)] <- 0
# xray_brick_sum <- as.data.frame(xray_brick, xy=T) %>%
#   select(x, y) %>%
#   bind_cols(xray_brick_summary)

directories <- list.dirs("/Users/Caitlin/Desktop/dataStitcher/example_dataset", full.names = T , recursive =F)
directories <- directories[!str_detect(directories, "Unknown|SEM|Os")]

overview_brick_list <- list()
overview_data <- list()

for(i in 1:length(directories)){
  path = directories[i]
  files <- list.files(path, full.names = T, pattern = "overview.*tif")
  if(length(files) >0){
    overview_data[[i]] <- str_extract(path, "([^/]+$)")
    message(paste0("Stacking ",overview_data[[i]], " data (element ", i, " of ", length(directories), ")..."))
    image <- files %>% image_read() %>% 
      image_quantize(colorspace = 'gray') %>% 
      image_equalize() 
    temp_file <- tempfile()
    image_write(image, path = temp_file, format = 'tiff')
    image <- raster(temp_file) %>%
      cut(breaks = c(-Inf, 150, Inf)) - 1
    image <- aggregate(image, fact=4)
    overview_brick_list[[i]] <- image
  }
}

message("Stacking complete. Creating x-ray brick...")
overview_brick <- do.call(brick, na.omit(overview_brick_list))
names(overview_brick) <- unlist(overview_data)
message("...complete.")


#bulk chemistry of overview area
overview_summary <- as.data.frame(overview_brick, xy = T)%>%
  select(-x, -y) %>%
  colSums(na.rm = T) %>%
  as.data.frame() %>%
  rename_(value = ".") %>%
  mutate(total = sum(value),
         value = (value/total)*100,
         element = names(overview_brick))


#bulk chemistry of transect 
xray_summary <- as.data.frame(xray_brick, xy = T)%>%
  select(-x, -y) %>%
  colSums(na.rm = T) %>%
  as.data.frame() %>%
  rename_(value = ".") %>%
  mutate(total = sum(value),
         value = (value/total)*100,
         element = names(xray_brick))


# spatial stats -----------------------------------------------------------
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

#create ppp for point pattern analysis in spatstat
window <- owin(xrange = c(st_bbox(cell_centroids)$xmin, st_bbox(cell_centroids)$xmax),
               yrange = c(st_bbox(cell_centroids)$ymin, st_bbox(cell_centroids)$ymax))
cell_centroids_ppp <- ppp(cell_centroids_coords$X, cell_centroids_coords$Y, window) %>%
  rescale(22.15, "μm")

#compute quadrat grid density and intensity -- add dynamic argument for nx/ny 

cell_dens_intensity <- as.data.frame(intensity(cells_quadrat, image=F), xy = T) %>%
  rename(intensity = Freq)

cells_quadrat <- as.data.frame(quadratcount(cell_centroids_ppp, nx= 12, ny=3), xy=T) %>%
  bind_cols(data.frame("grid_id" = unlist(rev(split(rev(1:36), rep_len(1:12, length(1:36))))))) %>%
  left_join(cell_dens_intensity) %>%
  select(Freq, grid_id, intensity) 

quadrat_grid <- st_make_grid(SEM_panorama_merged, cellsize = c(442, 248)) %>% 
  st_sf(grid_id = 1:length(.)) %>% left_join(cells_quadrat)

grid_lab <- st_centroid(quadrat_grid) %>% cbind(st_coordinates(.)) %>% left_join(cells_quadrat)

# view the sampled points, polygons and grid
rasterVis::gplot(SEM_panorama_merged) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'gray', high = 'white') +
  coord_fixed() +
  ggnewscale::new_scale_fill() +
  geom_sf(data = quadrat_grid, inherit.aes = F, aes(fill = intensity), alpha = 0.5, lwd = 0.3, color = "black") +
  scale_fill_viridis_c() +
  geom_sf(data = cells_polygon_projected, inherit.aes = F, fill = "red", lwd = 0) + 
  geom_text(data = grid_lab, aes(x = X, y = Y, label = Freq), size = 3, color = "black", fontface = "bold") +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") 
  

##point process model

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

#more resources 
#https://eburchfield.github.io/files/Point_pattern_LAB.html



