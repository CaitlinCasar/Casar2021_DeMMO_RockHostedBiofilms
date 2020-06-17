suppressPackageStartupMessages(require(optparse))

# opt parse arguments -----------------------------------------------------

option_list = list(
  make_option(c("-f", "--file"), action="store", default=getwd(), type='character',
              help="Directory containing files from dataStitchR output and feature data."),
  make_option(c("-k", "--transect"), action="store", default=NA, type='character',
              help="Regex pattern for brick of xray transect from dataStitchR output"),
  make_option(c("-t", "--overview"), action="store", default=NA, type='character',
              help="Optional regex pattern of brick of overview area from dataStitchR output"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="Output file directory. This is where your reports and plots will be saved."),
  make_option(c("-n", "--name"), action="store", default=NA, type='character',
              help="Optional name for output files."),
  make_option(c("-b", "--base"), action="store", default=NA, type='character',
              help="Required file path for SEM panoramic tif from dataStitchR ouput."),
  make_option(c("-j", "--overview_base"), action="store", default=NA, type='character',
              help="Optional file path for overview SEM tif."),
  make_option(c("-c", "--cores"), action="store", default=1, type='numeric',
              help="Number of cores to use for parallel processing, default is 1."),
  make_option(c("-u", "--use-scale"), action="store", default=1, type='numeric',
              help="Optional scale to use for converting pixels to microns."),
  make_option(c("-z", "--cells"), action="store", default=NA, type='character',
              help="Required file path to tif file of cell features."),
  make_option(c("-a", "--biogenic"), action="store", default=NA, type='character',
              help="Optional filt path for tif file of biogenic features."),
  make_option(c("-m", "--model_vars"), action="store", default=1, type='numeric',
              help="Required number of variables to include in point process models where variables are elements."),
  make_option(c("-d", "--d_cols"), action="store", default=12, type='numeric',
              help="Optional number of columns for quadrat grid."),
  make_option(c("-y", "--y_rows"), action="store", default=3, type='numeric',
              help="Optional number of rows for quadrat grid."),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print updates to console [default %default]."),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print anything to the console."),
  make_option(c("-p", "--pdf"), action="store", default=FALSE,
              help="Generate element plot with density contours, default is FALSE [default %default].")  
)
opt = parse_args(OptionParser(option_list=option_list))



#load dependencies
pacman::p_load(MASS, parallel, spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, units, raster, rgeos, imager, ggnewscale,  cowplot, tidyverse, rgdal, rasterVis)
#removed fasterRaster, lwgeom,stars,

# import the data ---------------------------------------------------------
message(paste0("Importing data from ", opt$n, "..."))
files <- list.files(opt$f, full.names = T)

#load xray raster bricks
message("importing xray brick")
xray_brick <- brick(files[str_detect(files, opt$transect)]) 

if(!is.na(opt$overview)){
  message("importing overview brick")
  overview_brick <- brick(files[str_detect(files, opt$overview)])
}
if(!is.na(opt$overview_base)){
  message("importing overview base SEM image")
  overview_SEM_image <- raster(files[str_detect(files, opt$overview_base)])
}
message("importing SEM pano image")
#load base SEM image
SEM_image <- raster(files[str_detect(files, opt$base)])

#cell and biogenic feature image paths
cells <- files[str_detect(files, opt$cells)]
if(!is.na(opt$biogenic)){
  biogenic <- files[str_detect(files, opt$biogenic)]
}

#set scale in number of pixels per micron 
micron_scale <- opt$u
#micron_scale <- 4.2

#create output directory
if(!is.na(opt$out)){
  dir.create(opt$out)
}

#plotting function
plot_PDF <- function(data_to_plot, plot_name){
  if(!is.na(opt$n)){
    if(!is.na(opt$out)){
      pdf(paste0(opt$out, "/", opt$n, "_", plot_name,".pdf"),
          width = 13.33, 
          height = 7.5)
    }else{
      pdf(paste0(opt$n, "_", plot_name,".pdf"),
          width = 13.33, 
          height = 7.5)
    }
  }else{
    if(!is.na(opt$out)){
      pdf(paste0(opt$out, "/", plot_name,".pdf"),
          width = 13.33, 
          height = 7.5)
    }else{
      pdf(paste0(plot_name,".pdf"),
          width = 13.33, 
          height = 7.5) 
    }
  }
  
  
  print(data_to_plot)
  
  dev.off()
}

#data writing function
write_data <- function(data, file_name){
  if(!is.na(opt$n)){
    if(!is.na(opt$out)){
      write_csv(data, paste0(opt$out, "/", opt$n, "_", file_name,".csv"))
    }else{
      write_csv(data, paste0(opt$n, "_", file_name,".csv"))
    }
  }else{
    if(!is.na(opt$out)){
      write_csv(data, paste0(opt$out, "/", file_name,".csv"))
    }else{
      write_csv(data, paste0(file_name,".csv"))
    }
  }
}


# rasterize/polygonize cells and biogenic features ------------------------

#function for rasterizing or polygonizing
image_to_polygon <- function(image_path, to_raster = TRUE, pres_abs = TRUE, drop_and_fill = TRUE, polygonize = TRUE){
    message("reading image...")
    image <- image_path 
    if(to_raster == T){
      message("rasterizing...")
      image_raster_noThresh <- raster(image) 
      image_raster <- image_raster_noThresh %>% setValues(if_else(values(image_raster_noThresh) > 0, 1, 0))
      image_raster 
      if(polygonize == TRUE){
        message("polygonizing...")
        image_polygon <- rasterToPolygons(image_raster, function(x){x == 1}, dissolve = TRUE) %>% 
          st_as_sf() 
        if(drop_and_fill == TRUE){
          message("processing drop and fill...")
          image_polygon_drop_fill <- image_polygon %>% 
            st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
            fill_holes(threshold = 50) %>%
            drop_crumbs(threshold = 10) %>% #raster to polygon was creating some tiny polygons 
            st_cast("POLYGON")
          image_polygon_drop_fill
        }else{
          image_polygon %>%
            st_set_crs("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m") %>%
            st_cast("POLYGON")
        }
      }else{
        image_raster
      }
    }else{
      image
    }
}

message("transforming cell data...")
#generate cell and biogenic feature polygons and rasters
cells_polygon <- image_to_polygon(cells, drop_and_fill = F)
cells_raster <- image_to_polygon(cells, polygonize = F)
if(!is.na(opt$biogenic) & length(biogenic) > 0){
  message("transforming biogenic data...")
  biogenic_polygon <- image_to_polygon(biogenic, drop_and_fill = F)
  biogenic_raster <- image_to_polygon(biogenic, polygonize = F)
}


# plot check for testing --------------------------------------------------
#check plotting if cells are polygons
# rasterVis::gplot(SEM_image) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = 'black', high = 'white') +
#   coord_fixed() +
#   geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
#   geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
#   coord_sf(datum = NA)

# calculate cell centroids for modeling  -----------------------------------------------------------
message("Preparing data for modeling...")
# calculate cell centroids 
cell_centroids <- st_centroid(cells_polygon)

cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame() 
#create ppp for point pattern analysis in spatstat
window <- owin(xrange = c(st_bbox(xray_brick)$xmin, st_bbox(xray_brick)$xmax),
               yrange = c(st_bbox(xray_brick)$ymin, st_bbox(xray_brick)$ymax))

cell_centroids_ppp <- ppp(cell_centroids_coords$X, cell_centroids_coords$Y, window) %>%
  rescale(micron_scale, "μm")


# plot check for testing --------------------------------------------------


#check plotting with density contours
# rasterVis::gplot(SEM_image) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = 'black', high = 'white') +
#   coord_fixed() +
#   ggnewscale::new_scale_color() +
#   geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
#   scale_color_viridis_c() +
#   geom_sf(data = cells_polygon, inherit.aes = F, fill = "orange", lwd = 0) +
#   geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
#   coord_sf(datum = NA) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# plot spatial randomness ppm models --------------------------------------------------

#do we expect to see clustering?
#Where K falls under the theoretical Kpois line the points are more clustered at distance r, and vis versa
message("Modeling spatial randomness...")
cell_fest <- Fest(cell_centroids_ppp) %>%
  dplyr::select(-hazard, -theohaz) %>%
  gather(type, value, theo:km) %>%
  mutate(test = "F function")

cell_kest <- Kest(cell_centroids_ppp) %>%
  gather(type, value, theo:iso) %>%
  mutate(test = "K function")

cell_jest <- Jest(cell_centroids_ppp) %>%
  dplyr::select(-hazard) %>%
  gather(type, value, theo:km) %>%
  mutate(test = "J function")

cell_gest <- Gest(cell_centroids_ppp) %>%
  dplyr::select(-hazard, -theohaz) %>%
  gather(type, value, theo:km) %>%
  mutate(test = "G function")

ppm_all <- cell_fest %>%
  bind_rows(cell_kest) %>%
  bind_rows(cell_jest) %>%
  bind_rows(cell_gest) 

ppm_all_plot <- ppm_all %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = value, color = type)) +
  xlab(expression("Distance ("~mu~"m)")) +
  facet_wrap(~test, scales = "free")


# generate PDF 

plot_PDF(ppm_all_plot, "ppm_models")

#write data
write_data(ppm_all, "ppm_models")


# Calculate ANN -----------------------------------------------------------
message("Calculating ANN...")
ANN <- apply(nndist(cell_centroids_ppp, k=1:100),2,FUN=mean)

ANN_plot <- ggplot() +
  geom_point(aes(x=1:100, y=ANN)) +
  xlab ("Nth nearest neighbor") +
  ylab(expression("Average distance ("~mu~"m)"))

plot_PDF(ANN_plot, "ANN")

ann.p <- mean(nndist(cell_centroids_ppp, k=1))

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

#plot random vs true points with density contours
ANN_dens_plot <- ggplot()+
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

plot_PDF(ANN_dens_plot, "ANN_density")

#plot ANN random vs. population
ANN_hist <- ggplot(as.data.frame(ann.r), aes(x=ann.r)) + 
  geom_histogram(color="#999999", bins=100, alpha = 0.5) +
  geom_vline(aes(xintercept=mean(ann.r)),
             color="#2b2d2f", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=ann.p),
             color="#E69F00", linetype="dashed", size=1) +
  xlim(0, max(ann.r) + 10) +
  labs(y= "Frequency", x = expression("ANN Distance ("~mu~"m)")) +
  annotate("text", x = 0, y =5, size = 3, label = paste0("pseudo p-value: ", "\n", p), hjust = "left")

plot_PDF(ANN_hist, "ANN_hist")


#p = chance that given an infinite number of simulations at least one realization of a point pattern could produce an ANN value more extreme than this sample

# element point process models --------------------------------------------
message("Modeling cell distribution as a function of rock chemistry...")
#test whether elements explain cell distribution

#add a limit to number of model variables based on rel abundance of elements

#generate list of all possible element combinations
element_combos = list()
for(element in 1:opt$model_vars){
  element_combos = append(element_combos, combn(names(xray_brick), element, simplify = F))
}

#generate images for each element
for(k in 1:length(names(xray_brick))){
    id = names(xray_brick)[k]
    assign(paste0(id), as.im(xray_brick[[id]]) %>% rescale(micron_scale, "μm"))
}

#calculate p-values for each element combination to see if any significantly expain cell distribution

element_models <- list()
for(i in 1:length(element_combos)){
  element_models[[i]] <- as.formula(paste("ppm(cell_centroids_ppp ~", paste(unlist(element_combos[i]), collapse = "+"), ")"))
}

#initialize empty lists for storing model outputs
element_covariate_models <- list()
element_model_results <- data.frame(model_ID = as.numeric(), model = as.character(), Npar=as.numeric(), name = as.character(), value = as.numeric()) 

#create the null model for comparison
ppm0 <- ppm(cell_centroids_ppp ~ 1) 

for(i in 1:length(element_models)){
  ppm1 <- ppm(element_models[[i]]) #create a model with covariates
  anova_lrt <- tidy(anova(ppm0,ppm1, test="LRT")) %>% 
    pivot_longer(cols = df:p.value, values_to = "value") %>%
    mutate(model_ID = i, 
           model = format(element_models[[i]]),
           n_elements = str_count(model, "[+]") + 1)
  element_model_results <- element_model_results %>%
    bind_rows(anova_lrt)
  element_covariate_models <- append(element_covariate_models, list(ppm1)) #store each model in the element_covariate_models list
}

#filter for only significant models based on anova results 
significant_models <- element_model_results %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = name, values_from=value) %>%
  filter(p.value <= 0.05)

#subset the full list of models for only the significant ones 
#sig_model_data <- element_covariate_models[significant_models$model_ID]


best_model <- significant_models %>%
  filter(Deviance == max(Deviance)) %>%
  filter(p.value == min(p.value)) %>%
  filter(n_elements == min(n_elements))

best_model_data <- element_covariate_models[best_model$model_ID]

saveRDS(best_model_data, file = paste0(opt$n,"_best_models.Rds"))

write_data(significant_models, "significant_models")

# plot Spearman local correlation over cells + elements from most --------
message("Modeling local correlation between cells and rock elements...")
sig_elements <- which( names(xray_brick) %in% (element_probs %>% filter(pval == min(na.omit(element_probs$pval))) %>% dplyr::select(elements) %>% unlist()))

element_raster_list <- list()
for(i in 1:length(sig_elements)){
  element_raster_list[i] <- xray_brick[[sig_elements[i]]]
}
cells_aggregated <- aggregate(cells_raster, 3, mean)

cell_element_corr_fun <- function(element){
  element_aggregated <- aggregate(element, 3, mean)
  corLocal(element_aggregated, cells_aggregated, method='spearman')
}
cell_element_corr <- lapply(element_raster_list, cell_element_corr_fun)

#stack rasters
cell_element_corr_stack <- stack(cell_element_corr)
names(cell_element_corr_stack) <- names(xray_brick)[sig_elements]

#create custom themes for levelplot
zeroCol <- NA
myTheme <- rasterTheme(region = c(zeroCol, viridis::viridis(30)))
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

#generate PDF
local_corr_plot <- background_image +  element_background_plot + cell_element_corr_plot
plot_PDF(local_corr_plot, "local_correlation")

# calculate total correlation between cells/biogenic features and elements--------
message("Calculating total correlation between cells/biogenic features and elements...")
if(!is.na(opt$biogenic) & length(biogenic) > 0){
  total_cor <- data.frame(element = character(), cell_cor = numeric(), cell_pval = numeric(), biogenic_cor = numeric(), biogenic_pval = numeric())
  cell_mat <- as.matrix(cells_raster)
  biogenic_mat <- as.matrix(biogenic_raster)
  for(i in 1:length(names(xray_brick))){
    element_mat <- as.matrix(xray_brick[[i]])
    result <- cor.test(element_mat, cell_mat)
    result2 <- cor.test(element_mat, biogenic_mat)
    total_cor <- total_cor %>% add_row(element = names(xray_brick[[i]]), 
                                       cell_cor = result$estimate, cell_pval = result$p.value,
                                       biogenic_cor = result2$estimate, biogenic_pval = result2$p.value)
    
  }
  #write data to csv
  write_data(total_cor, "total_correlation")
  
  #generate PDF
  total_cor <- total_cor %>%
    dplyr::select(-cell_pval, -biogenic_pval) %>%
    gather(type, cor, cell_cor:biogenic_cor) %>%
    ggplot() +
    geom_line(aes(reorder(element, cor), cor, color = type, group=type)) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    xlab("Element") +
    ylab("Correlation")
  
  plot_PDF(total_cor, "total_correlation")
}else{
  total_cor <- data.frame(element = character(), cell_cor = numeric(), cell_pval = numeric())
  cell_mat <- as.matrix(cells_raster)
  for(i in 1:length(names(xray_brick))){
    element_mat <- as.matrix(xray_brick[[i]])
    result <- cor.test(element_mat, cell_mat)
    total_cor <- total_cor %>% add_row(element = names(xray_brick[[i]]), 
                                       cell_cor = result$estimate, cell_pval = result$p.value)
    
  }
  #write data to csv
  write_data(total_cor, "total_correlation")
  
  #generate PDF
  total_cor <- total_cor %>%
    ggplot() +
    geom_line(aes(reorder(element, cell_cor), cell_cor, group=1)) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    xlab("Element") +
    ylab("Correlation")
  
  plot_PDF(total_cor, "total_correlation")
}





# generate cell stats -----------------------------------------------------
message("Generating cell summary report...")
#function to draw ellipses around cell polygons, measure major and minor axes 
ellipsoid_fun <- function(data){
  polygon_coords <- data %>%
    st_coordinates() %>%
    as.data.frame() %>%
    dplyr::select(-L1, -L2) %>%
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
  mutate(area = st_area(geometry)/(micron_scale^2))#,
         #perimeter = st_perimeter(geometry)/micron_scale)
cells_polygon_stats$roundness <- sapply(cells_polygon$geometry, ellipsoid_fun)
cells_polygon_stats <- cells_polygon_stats %>%
  mutate(shape = case_when(roundness < 1.12 ~ "cocci",
                           roundness >= 1.12 & roundness < 5 ~"rod",
                           TRUE ~ "filament"))
if(!is.na(opt$biogenic) & length(biogenic) > 0){
  biogenic_polygon_stats <- biogenic_polygon %>%
    mutate(area = st_area(geometry)/(micron_scale^2))
}
  

transect_area <- ((extent(SEM_image)[2])/micron_scale)*((extent(SEM_image)[4])/micron_scale)

cell_summary <- cells_polygon_stats %>% 
  as.data.frame %>%
  group_by(shape) %>% 
  summarise(value = n()) %>%
  rename(observation = shape) %>%
  bind_rows(data.frame(observation = "mean cell area (μm^2)", value = mean(cells_polygon_stats$area)),
            data.frame(observation = "total cell area (μm^2)", value = sum(cells_polygon_stats$area)),
            data.frame(observation = "coverage cell area (%)", value = ((sum(cells_polygon_stats$area))/transect_area)*100),
            data.frame(observation = "total cells", value = nrow(cells_polygon_stats)),
            data.frame(observation = "cell density (cells/cm^2)", value = nrow(cells_polygon_stats)/(transect_area/100000000)),
            data.frame(observation = "cell ANN", value = ann.p),
            data.frame(observation = "mean random ANN", value = mean(ann.r)),
            data.frame(observation = "pseudo p-value", value=p))
            
if(!is.na(opt$biogenic) & length(biogenic) > 0){
  cell_summary <- cell_summary %>%
  bind_rows(data.frame(observation = "total biogenic area (μm^2)", value = sum(biogenic_polygon_stats$area)),
  data.frame(observation = "coverage biogenic area (%)", value = ((sum(biogenic_polygon_stats$area))/transect_area)*100))
}            

#write data to csv
write_data(cell_summary, "cell_summary")

# #generate chemistry stats -----------------------------------------------
message("Generating rock chemistry report...")
#bulk chemistry of transect 
transect_summary <- as.data.frame(xray_brick, xy = T) %>%
  replace(is.na(.), 0) %>%
  summarise_all(funs(sum)) %>%
  gather(element, `wt%`, -x, -y) %>%
  dplyr::select(-x, -y) %>% 
  mutate(`wt%` = `wt%`/sum(`wt%`)*100,
         type = "transect") %>%
  spread(type, `wt%`)

#bulk chemistry of overview area
if(!is.na(opt$overview)){
  overview_summary <- as.data.frame(overview_brick, xy = T) %>%
    replace(is.na(.), 0) %>%
    summarise_all(funs(sum)) %>%
    gather(element, `wt%`, -x, -y) %>%
    dplyr::select(-x, -y) %>% 
    mutate(`wt%` = `wt%`/sum(`wt%`)*100,
           type = "overview") %>%
    spread(type, `wt%`)
  element_summary <- overview_summary %>%
    full_join(transect_summary)
}else{
  element_summary <- transect_summary
}

#write data to csv
write_data(element_summary, "element_summary")
# plot element data with density contours -----------------------------------------------------------

#generate element map with cells 
if(opt$p){
message("Generating density contour plot...")
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
    p + 
      coord_fixed() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none")
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
#element_plot <- element_plotter(xray_frame, xray_brick, SEM_image, element_colors, density = F)

##element plot with density contours
element_plot <- element_plotter(xray_frame, xray_brick, SEM_image, element_colors)
SEM_plot <- rasterVis::gplot(SEM_image) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  coord_fixed() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
density_plot <- rasterVis::gplot(SEM_image) +
  geom_tile(aes(fill = value), alpha = 0.2) +
  scale_fill_gradient(low = 'black', high = 'white') +
  ggnewscale::new_scale_fill()

if(!is.na(opt$biogenic) & length(biogenic) > 0){
  density_plot <- density_plot + 
    coord_fixed() +
    ggnewscale::new_scale_color() +
    geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
    scale_color_viridis_c() +
    geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
    geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0) +
    coord_sf(datum = NA)  +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
}else{
  density_plot <- density_plot + 
    coord_fixed() +
    ggnewscale::new_scale_color() +
    geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
    scale_color_viridis_c() +
    geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
    coord_sf(datum = NA)  +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
}

element_plot_with_legend <- plot_grid(
  SEM_plot,
  density_plot,
  element_plot, 
  plot_grid(get_legend(element_plot_legend), 
            ncol = 1), 
  nrow = 4, 
  rel_heights = c(4, 4, 4, 2)
)


#generatePDF
plot_PDF(element_plot_with_legend, "density_contour_with_elements")

}

# element overview plot ---------------------------------------------------

if(!is.na(opt$overview) & !is.na(opt$overview_base)){
  message("geerating overview element plot...")
  overview_xray_frame <- as.data.frame(overview_brick, xy=TRUE) 
  overview_xray_frame <- gather(overview_xray_frame, element, value, colnames(overview_xray_frame)[3]:colnames(overview_xray_frame)[ncol(overview_xray_frame)])

  element_plot <- element_plotter(overview_xray_frame, overview_brick, overview_SEM_image, element_colors)
  SEM_plot <- rasterVis::gplot(overview_SEM_image) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'black', high = 'white') +
    coord_fixed() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
  element_plot_legend <- data.frame(element = unique(overview_xray_frame$element)) %>%
    rownames_to_column() %>% 
    ggplot(aes(element, rowname, fill=element)) + 
    geom_bar(stat= "identity") + 
    scale_fill_manual(values = element_colors) +
    guides(fill=guide_legend(nrow=2)) +
    theme(legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 8))
  
  element_plot_with_legend <- plot_grid(
    SEM_plot,
    element_plot, 
    plot_grid(get_legend(element_plot_legend), 
              ncol = 1), 
    nrow = 3, 
    rel_heights = c(4, 4, 2)
  )
  
  
  #generatePDF
  plot_PDF(element_plot_with_legend, "overview_with_elements")
  
}
# quadrat plot ------------------------------------------------------------
message("Generating quadrat plot...")

#compute quadrat grid density and intensity  
quad_nx = opt$d
quad_ny = opt$y

cells_quadrat <- quadratcount(cell_centroids_ppp, nx = quad_nx, ny = quad_ny)
cell_dens_intensity <- as.data.frame(intensity(cells_quadrat, image=F), xy = T) %>%
  rename(intensity = Freq)

cells_quadrat <- as.data.frame(quadratcount(cell_centroids_ppp, nx = quad_nx , ny = quad_ny), xy=T) %>%
  bind_cols(data.frame("grid_id" = unlist(rev(split(rev(1:(quad_nx*quad_ny)), rep_len(1:quad_nx, length(1:(quad_nx*quad_ny)))))))) %>%
  left_join(cell_dens_intensity) %>%
  dplyr::select(Freq, grid_id, intensity) 

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
  scale_fill_viridis_c() +
  geom_sf(data = cells_polygon, inherit.aes = F, fill = "red", lwd = 0) + 
  geom_text(data = grid_lab, aes(x = X, y = Y, label = Freq), size = 3, color = "black", fontface = "bold") +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") 

#generate PDF
plot_PDF(quadrat_plot, "quadrat_plot")

message("...complete.")
