suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-f", "--file"), action="store", default=getwd(), type='character',
              help="Name of brick file."),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="Output file directory. This is where your x-ray raster brick and output figures will be saved."),
  make_option(c("-n", "--name"), action="store", default=NA, type='character',
              help="Optional name for output files."),
  make_option(c("-b", "--base-images"), action="store", default=NA, type='character',
              help="SEM image file directory."),
  make_option(c("-c", "--coords"), action="store", default=NA, type='character',
              help="Tab-delimited file of xy coordinates for each image. A third column should denote stitching positions that correspond to the file names for each image."),
  make_option(c("-u", "--use-positions"), action="store", default="-?(?<![Kα1||Kα1_2])\\d+", type='character',
              help="Optional regex pattern to extract position IDs from each file name that corresponds to positions in the xy file. The default searches for numbers that appear after 'Kα1' or 'Kα2'. Numbers can include signs, i.e. -1 is acceptable."),
  make_option(c("-z", "--z-format"), action="store", default="*", type='character',
              help="Optional regex pattern of x-ray image formats to select for stitching, i.e. '.tif'."),
  make_option(c("-m", "--make"), action="store", default="*", type='character',
              help="Optional regex pattern of SEM image formats to select for stitching, i.e. '.tif'. You do not need to specify this unless you are generating a PDF output."),
  make_option(c("-a", "--all-exclude"), action="store", default=NA, type='character',
              help="Optional regex pattern of x-ray file directories to exclude from stitiching, i.e. the element your sample was coated with."),
  make_option(c("-d", "--drop"), action="store", default=NA, type='character',
              help="Optional regex pattern of files to exclude from x-ray data stitching."),
  make_option(c("-y", "--y-exclude"), action="store", default=NA, type='character',
              help="Optional regex pattern of files to exclude from SEM image stitiching. You do not need to specify this unless you are generating a PDF output."),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print updates to console [default %default]."),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print anything to the console."),
  make_option(c("-p", "--pdf"), action="store", default=FALSE,
              help="Generate PDF of x-ray brick colored by element superimposed on the SEM image, default is TRUE [default %default].")  
)
opt = parse_args(OptionParser(option_list=option_list))



#load dependencies
pacman::p_load(spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, lwgeom, units, raster, rgeos, imager,ggnewscale,  magick, stars, fasterRaster, cowplot, tidyverse, rgdal, rasterVis)


#set working directory
setwd("/Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO3/D3T13exp_Dec2019_Poorman/")

#load xray raster bricks
xray_brick <- brick("D3T13exp_Dec2019_brick.grd") 
overview_brick <- brick("D3T13exp_Dec2019_overview_brick.grd")

#load base SEM image
SEM_image <- raster("D3T13exp_Dec2019_SEM_pano.tif")

#cell and biogenic feature image paths
cells <- "cells.tif"
biogenic <- "biogenic.tif"

#set scale in number of pixels per micron 
micron_scale <- 4.2

#function for rasterizing or polygonizing
image_to_polygon <- function(image_path, to_raster = TRUE, pres_abs = TRUE, drop_and_fill = TRUE, polygonize = TRUE, equalizer = TRUE){
  if(equalizer == T){
    message("reading image...")
    image <- image_path %>% image_read() %>% image_quantize(colorspace = 'gray') %>% image_equalize() 
    image
    if(to_raster == T){
      message("rasterizing...")
      temp_file <- tempfile()
      image_write(image, path = temp_file, format = 'tiff')
      image_raster_noThresh <- raster(temp_file) #%>%
        #cut(breaks = c(-Inf, 150, Inf)) - 1
      image_raster <- image_raster_noThresh %>% setValues(if_else(values(image_raster_noThresh) > 150, 0, 1))
      image_raster 
      if(polygonize == TRUE){
        message("polygonizing...")
        image_polygon <- rasterToPolygons(image_raster, function(x){x == 1}, dissolve = TRUE) %>% 
          st_as_sf() 
        if(drop_and_fill == TRUE){
          message("processing drop and fill...")
          image_polygon_drop_fill <- image_polygon %>% 
            st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
            fill_holes(threshold = 100) %>%
            drop_crumbs(threshold = 50) %>% #raster to polygon was creating some tiny polygons 
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

#generate cell and biogenic feature polygons and rasters
cells_polygon <- image_to_polygon(cells)
biogenic_polygon <- image_to_polygon(biogenic)

cells_raster <- image_to_polygon(cells, polygonize = F)
biogenic_raster <- image_to_polygon(biogenic, polygonize = F)

#check plotting if cells are polygons
rasterVis::gplot(SEM_image) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  coord_fixed() +
  geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
  geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
  coord_sf(datum = NA)


# calculate cell centroids 
cell_centroids <- st_centroid(cells_polygon)

cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame() #%>% #this code comes before cells_polygon_stats defined
  #mutate(shape = cells_polygon_stats$shape,
         #area = cells_polygon_stats$area) 

#check plotting with density contours
rasterVis::gplot(SEM_image) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  coord_fixed() +
  ggnewscale::new_scale_color() +
  geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
  scale_color_viridis_c() +
  geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
  geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
  coord_sf(datum = NA) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")


# spatial stats -----------------------------------------------------------
message("Calculating statistics...")



# quadrat stats -----------------------------------------------------------

#create ppp for point pattern analysis in spatstat
window <- owin(xrange = c(st_bbox(cell_centroids)$xmin, st_bbox(cell_centroids)$xmax),
               yrange = c(st_bbox(cell_centroids)$ymin, st_bbox(cell_centroids)$ymax))

cell_centroids_ppp <- ppp(cell_centroids_coords$X, cell_centroids_coords$Y, window) %>%
  rescale(micron_scale, "μm")

#compute quadrat grid density and intensity -- add dynamic argument for nx/ny 
quad_nx = 12
quad_ny = 3

cells_quadrat <- quadratcount(cell_centroids_ppp, nx = quad_nx, ny = quad_ny)
cell_dens_intensity <- as.data.frame(intensity(cells_quadrat, image=F), xy = T) %>%
  rename(intensity = Freq)

#*******make nx ny dynamic*********

cells_quadrat <- as.data.frame(quadratcount(cell_centroids_ppp, nx = quad_nx , ny = quad_ny), xy=T) %>%
  bind_cols(data.frame("grid_id" = unlist(rev(split(rev(1:(quad_nx*quad_ny)), rep_len(1:quad_nx, length(1:(quad_nx*quad_ny)))))))) %>%
  left_join(cell_dens_intensity) %>%
  select(Freq, grid_id, intensity) 
  


# point process modeling --------------------------------------------------

#calculate spatial randomness 
#do we expect to see clustering?
#Where K falls under the theoretical Kpois line the points are more clustered at distance r, and vis versa

cell_kest <- Kest(cell_centroids_ppp)

kest_plot <- cell_kest %>%
  gather(type, value, theo:iso) %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = value, color = type))

# generate PDF ------------------------------------------------------------
message("Generating Kest plot...")

sample_id <- "D1T1exp_Dec2019_Poorman"
pdf(paste0(sample_id, "_Kest_plot.pdf"),
    width = 13.33, 
    height = 7.5)

print(element_plot_with_legend)

dev.off()

message("...complete.")


# Calculate ANN -----------------------------------------------------------

ANN <- apply(nndist(cell_centroids_ppp, k=1:100),2,FUN=mean)
plot(ANN ~ eval(1:100), type="b", main=NULL, las=1)


ann.p <- mean(nndist(cell_centroids_ppp, k=1))

#element_probabilities <- list()  ##not sure what I was going for here...
#generate random distribution of cells as null model
n = 599L
ann.r <- vector(length=n) # Create an empty object to be used to store simulated ANN values
for (i in 1:599L){
  rand.p   <- rpoint(n=cell_centroids_ppp$n, win=window) 
  rand.p <- ppp(rand.p$x, rand.p$y, window) %>%
    rescale(micron_scale, "μm")
  # Generate random point locations
  ann.r[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}

N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)


ggplot()+
  geom_point(aes(x = rand.p$x, y = rand.p$y), alpha=0.2, stroke=0) +
  geom_point(aes(x=cell_centroids_ppp$x, y=cell_centroids_ppp$y), color= "#E69F00", stroke=0) +
  coord_fixed()+
  ggnewscale::new_scale_color() +
  geom_density_2d(inherit.aes = F, mapping = aes(cell_centroids_ppp$x, cell_centroids_ppp$y, col = stat(level)/max(stat(level)))) +
  scale_color_viridis_c() +
  ggnewscale::new_scale_color() +
  geom_density_2d(inherit.aes = F, mapping = aes(rand.p$x, rand.p$y, col = stat(level)/max(stat(level))), alpha = 0.2) +
  scale_color_gradient(high = "black", low = "white") +
  #scale_color_grey(start = 0.8, end = 0.2)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#plot ANN random vs. population
ggplot(as.data.frame(ann.r), aes(x=ann.r)) + 
  geom_histogram(color="#999999", bins=100, alpha = 0.5) +
  geom_vline(aes(xintercept=mean(ann.r)),
             color="#2b2d2f", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=ann.p),
             color="#E69F00", linetype="dashed", size=1) +
  xlim(0, max(ann.r) + 10) +
  labs(y= "Frequency", x = "ANN Distance (μm)") +
  annotate("text", x = 0, y =5, size = 3, label = paste0("pseudo p-value: ", "\n", p), hjust = "left")


#p = chance that given an infinite number of simulations at least one realization of a point pattern could produce an ANN value more extreme than this sample


#element_probabilities[[length(names(xray_brick)) + 1]] <- min(N.greater + 1, n + 1 - N.greater) / (n +1)


# #test whether element distribution explains cell distribution 
# element_probabilities <- list()
# #this loop takes ~6 minutes per element...
# for(j in 1:length(names(xray_brick))){
#   message(paste0("Calculating covariate stats for ", names(xray_brick)[j], ", element ", j, " of ", length(names(xray_brick)), "..."))
#   n = 599L
#   ann.r <- vector(length=n)
#   xray_im <- as.im(xray_brick[[j]]) %>%
#     rescale(micron_scale, "μm")
#   for (i in 1:n){
#     rand.p   <- rpoint(n=cell_centroids_ppp$n, f = xray_im) 
#     ann.r[i] <- mean(nndist(rand.p, k=1))
#   }
#   N.greater <- sum(ann.r > ann.p)
#   element_probabilities[[j]] <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
# }

#point process model
#test whether elements explain cell distribution

#add a limit to number of model variables based on rel abundance of elements
n_variables <- nrow(element_summary %>% filter(transect >= 10)) #set threshold to 10%

#generate list of all possible element combinations
element_combos = list()
for(element in 1:n_variables){
  element_combos = append(element_combos, combn(names(xray_brick), element, simplify = F))
}

#generate images for each element
for(k in 1:length(names(xray_brick))){
    id = names(xray_brick)[k]
    assign(paste0(id), as.im(xray_brick[[id]]) %>% rescale(micron_scale, "μm"))
}

#calculate p-values for each element combination to see if any significantly expain cell distribution
# element_probabilities <- list()
# for(i in 1:length(element_combos)){
#   #print(paste(unlist(element_combos[i]), collapse = "+"))
#   ppm1 <- eval(parse(text=(paste("ppm(cell_centroids_ppp ~", paste(unlist(element_combos[i]), collapse = "+"), ")"))))
#   ppm0 <- ppm(cell_centroids_ppp ~ 1)
#   element_probabilities[[i]] <- unlist(anova(ppm1,ppm0, test="LRT")[4])[2]
# }

#this took ~20 minutes to run ~3000 out of 16k combos

element_models <- list()
for(i in 1:length(element_combos)){
  element_models[[i]] <- paste("ppm(cell_centroids_ppp ~", paste(unlist(element_combos[i]), collapse = "+"), ")")
}

element_prob_function <- function(element_model_list){
  ppm1 <- eval(parse(text=(element_model_list)))
  ppm0 <- ppm(cell_centroids_ppp ~ 1)
  unlist(anova(ppm1,ppm0, test="LRT")[4])[2]
}

#calculate probabilities in parallel using max # of detected cores 
library(parallel)
library(MASS) #note that this masks 'select' from tidyverse, 'area'amd 'select' from raster, amd 'area' from spatstat. Should be loaded before these packages or call these explicitely  
element_probabilities <- mclapply(element_models, element_prob_function, mc.cores = detectCores())


#my computer has 4 cores, so this will prob still take an hour to run 

element_probs <- data.frame("pval" = unlist(element_probabilities)) %>%
  mutate(model = element_models,
         elements = element_combos) 

sig_elements <- which( names(xray_brick) %in% (element_probs %>% filter(pval == min(na.omit(element_probs$pval))) %>% dplyr::select(elements) %>% unlist()))


#calculate Spearman correlation over pixels between cells + all elements and between biogenic + all elements
r1 <- aggregate(xray_brick[[16]], 3, mean)
r2 <- aggregate(cells_raster, 3, mean)
rc <- corLocal(r1, r2, method='spearman') #this step takes ~1 minute to run at agg factor of 3


element_raster_list <- list()
for(i in 1:length(sig_elements)){
  element_raster_list[i] <- xray_brick[[sig_elements[i]]]
}


cells_aggregated <- aggregate(cells_raster, 3, mean)

cell_element_corr_fun <- function(element){
  #message(paste0("Calculating Spearman correlation between cells and ", names(element), "..."))
  element_aggregated <- aggregate(element, 3, mean)
  corLocal(element_aggregated, cells_aggregated, method='spearman')
}
cell_element_corr <- lapply(element_raster_list, cell_element_corr_fun)
cell_element_corr_stack <- stack(cell_element_corr)
names(cell_element_corr_stack) <- element_probs %>% filter(pval == min(na.omit(element_probs$pval))) %>% dplyr::select(elements) %>% unlist()



#parallelizing does not work on this for some reason
#cell_element_corr <- mclapply(element_raster_list, cell_element_corr_fun, mc.cores = detectCores())

myTheme <- rasterTheme(region = c(zeroCol, viridis::viridis(30)))
#myTheme$panel.background$col = 'gray' 
zeroCol <- NA
element_theme <- rasterTheme(region = c(zeroCol, 'orange'))

# Customize the colorkey
my.at <- seq(-1, 1, length.out=length(myTheme$regions$col)-1)
my.ckey <- list(at=my.at, col=myTheme$regions$col)


background_stack <- list()
for(i in 1:length(sig_elements)){
  background_stack[[i]] <- SEM_image
}
background_stack <- stack(background_stack)
names(background_stack) <- names(cell_element_corr_stack)
background_image <- levelplot(background_stack, par.settings = GrTheme,
                              xlab=NULL, ylab=NULL, scales=list(draw=FALSE), alpha.regions = 0.5, colorkey=my.ckey)
background_image <- update(background_image, aspect=nrow(SEM_image)/ncol(SEM_image)) 
element_background_plot <- rasterVis::levelplot(reclassify(xray_brick[[sig_elements]], cbind(0, NA)), par.settings = element_theme,
                                               xlab=NULL, ylab=NULL, scales=list(draw=FALSE)) 
element_background_plot <- update(element_background_plot, aspect=nrow(SEM_image)/ncol(SEM_image)) 
cell_element_corr_plot <- rasterVis::levelplot(cell_element_corr_stack, par.settings = myTheme,
                     xlab=NULL, ylab=NULL, scales=list(draw=FALSE)) 
cell_element_corr_plot <- update(cell_element_corr_plot, aspect=nrow(SEM_image)/ncol(SEM_image)) #set the aspect ratio for the plots

background_image +  element_background_plot + cell_element_corr_plot

# test <- mask(xray_brick[[16]], cells_polygon)
# 
# 
# cor.test(as.matrix(test), as.matrix(cells_raster))
# cor_stack <- stack(xray_brick[[16]], cells_raster)
# x <- as.matrix(cor_stack)
# cor(x, method="spearman")
# 
# cell_element_cor <- list()
# for(i in 1:length(names(xray_brick))){
#   r1 <- as.matrix(xray_brick[[i]])
#   r2 <- as.matrix(cells_raster)
#   cell_element_cor[[i]] <- cor.test(r1, r2)
# }




#rescale raster layer 
range01 <- function(x){(na.omit(x)-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x)))}
test <- xray_brick[[1]] %>% setValues(range01(values(xray_brick[[1]])))

# summarize data ----------------------------------------------------------
message("Generating reports...")


# generate cell stats -----------------------------------------------------

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
  mutate(area = st_area(geometry)/(micron_scale^2),
         perimeter = st_perimeter(geometry)/micron_scale)
cells_polygon_stats$roundness <- sapply(cells_polygon$geometry, ellipsoid_fun)
cells_polygon_stats <- cells_polygon_stats %>%
  mutate(shape = case_when(roundness < 1.12 ~ "cocci",
                           roundness >= 1.12 & roundness < 5 ~"rod",
                           TRUE ~ "filament"))

cell_summary <- cells_polygon_stats %>% 
  as.data.frame %>%
  group_by(shape) %>% 
  summarise(value = n()) %>%
  rename(observation = shape) %>%
  bind_rows(data.frame(observation = "mean cell area (μm)", value = mean(cells_polygon_stats$area)),
            data.frame(observation = "total cells", value = nrow(cells_polygon_stats)),
            data.frame(observation = "cell density (cells/mm^2)", value = nrow(cells_polygon_stats)/((240*33.66516)/1000000)),
            data.frame(observation = "ANN", value = ann.p))



# #generate chemistry stats -----------------------------------------------

#bulk chemistry of overview area
overview_summary <- as.data.frame(overview_brick, xy = T) %>%
  replace(is.na(.), 0) %>%
  summarise_all(funs(sum)) %>%
  gather(element, `wt%`, -x, -y) %>%
  dplyr::select(-x, -y) %>% 
  mutate(`wt%` = `wt%`/sum(`wt%`)*100,
         type = "overview") %>%
  spread(type, `wt%`)

#bulk chemistry of transect 
transect_summary <- as.data.frame(xray_brick, xy = T) %>%
  replace(is.na(.), 0) %>%
  summarise_all(funs(sum)) %>%
  gather(element, `wt%`, -x, -y) %>%
  dplyr::select(-x, -y) %>% 
  mutate(`wt%` = `wt%`/sum(`wt%`)*100,
         type = "transect") %>%
  spread(type, `wt%`)

element_summary <- overview_summary %>%
  full_join(transect_summary)

# plot the data -----------------------------------------------------------

#generate element map with cells 

message("Generating element plot...")
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


element_plotter<-function(coord_frame, brick, SEM_image, colors, density=TRUE){
  p <-rasterVis::gplot(SEM_image) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'black', high = 'white') +
    ggnewscale::new_scale_fill()
  for(i in names(brick)){ 
    message(paste0("Adding ", names(brick[[i]]), " to plot..."))
    element_coords <- coord_frame %>%
      filter(element == names(brick[[i]]) & value!=0)
    p <- p+geom_raster(element_coords, mapping = aes(x, y, fill = element, alpha = value)) +
      scale_fill_manual(values = colors) +
      ggnewscale::new_scale_fill()
  }
  message("Writing element plot...")
  if(density){
    suppressWarnings(print(p + 
                             coord_fixed() +
                             ggnewscale::new_scale_color() +
                             geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
                             scale_color_viridis_c() +
                             geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
                             geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
                             coord_sf(datum = NA)  +
                             theme(axis.title = element_blank(),
                                   axis.text = element_blank(),
                                   legend.position = "none")))
  }else{
    suppressWarnings(print(p + 
                             coord_fixed() +
                             coord_sf(datum = NA)  +
                             geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
                             theme(axis.title = element_blank(),
                                   axis.text = element_blank(),
                                   legend.position = "none")))
  }
  
}

element_plot_legend <- data.frame(element = unique(xray_frame$element)) %>%
  rownames_to_column() %>% 
  ggplot(aes(element, rowname, fill=element)) + 
  geom_bar(stat= "identity") + 
  scale_fill_manual(values = element_colors) +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8))

##element plot with no density contours
element_plot <- element_plotter(xray_frame, xray_brick, SEM_image, element_colors, density = F)

##element plot with density contours
element_plot <- element_plotter(xray_frame, xray_brick, SEM_image, element_colors)

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

#generate quadrat plot
message("Generating quadrat plot...")

quadrat_grid <- st_make_grid(SEM_image, cellsize = c(extent(SEM_image)[2]/quad_nx, extent(SEM_image)[4]/quad_ny)) %>% 
  st_sf(grid_id = 1:length(.)) %>% left_join(cells_quadrat)

quadrat_grid_poly <- as.data.frame(st_coordinates(quadrat_grid)) %>%
  mutate(grid_id = rep(quadrat_grid$grid_id, 1, each=5)) %>%
  #filter(!row_number() %% 5 == 0) %>%
  left_join(quadrat_grid)

grid_lab <- st_centroid(quadrat_grid) %>% cbind(st_coordinates(.)) %>% left_join(cells_quadrat)

# view the sampled points, polygons and grid
quadrat_plot <- rasterVis::gplot(SEM_image) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'gray', high = 'white') +
  coord_fixed() +
  ggnewscale::new_scale_fill() +
  geom_polygon(data = quadrat_grid_poly, aes(X,Y, fill = intensity, group=grid_id), alpha = 0.5, lwd = 0.3, color = "black") +
  #geom_sf(data = quadrat_grid, inherit.aes = F, aes(fill = intensity), alpha = 0.5, lwd = 0.3, color = "black") + # this was causing issues with aspect ratio
  scale_fill_viridis_c() +
  geom_sf(data = cells_polygon, inherit.aes = F, fill = "red", lwd = 0) + 
  geom_text(data = grid_lab, aes(x = X, y = Y, label = Freq), size = 3, color = "black", fontface = "bold") +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") 

# generate PDF ------------------------------------------------------------
sample_id <- "D1T1exp_Dec2019_Poorman"
pdf(paste0(sample_id, "_quadrat_plot.pdf"),
    width = 13.33, 
    height = 7.5)

print(quadrat_plot)

dev.off()

message("...complete.")
