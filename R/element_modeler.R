
#load dependencies
pacman::p_load(MASS, parallel, spatstat, geostatsp, maptools, cluster, stringr, smoothr, sf, units, raster, rgeos, imager, ggnewscale,  cowplot, tidyverse, rgdal, rasterVis, broom)
#removed fasterRaster, lwgeom,stars,

# import the data ---------------------------------------------------------
message(paste0("Importing data from ", samplename, "..."))
files <- list.files(filename, full.names = T)

#load xray raster bricks
message("importing xray brick")
xray_brick <- brick(files[str_detect(files, transectname)]) 

if(!is.na(overviewname)){
  message("importing overview brick")
  overview_brick <- brick(files[str_detect(files, overviewname)])
}
if(!is.na(overview_base)){
  message("importing overview base SEM image")
  overview_SEM_image <- raster(files[str_detect(files, overview_base)])
}
message("importing SEM pano image")
#load base SEM image
SEM_image <- raster(files[str_detect(files, base_name)])

#cell and biogenic feature image paths
cells <- files[str_detect(files, cellsname)]
if(!is.na(biogenicname)){
  biogenic <- files[str_detect(files, biogenicname)]
}

#set scale in number of pixels per micron 
micron_scale <- scale_conversion
#micron_scale <- 4.2

#create output directory
if(!is.na(outname)){
  dir.create(outname)
}

#plotting function
plot_PDF <- function(data_to_plot, plot_name){
  if(!is.na(samplename)){
    if(!is.na(outname)){
      pdf(paste0(outname, "/", samplename, "_", plot_name,".pdf"),
          width = 13.33, 
          height = 7.5)
    }else{
      pdf(paste0(samplename, "_", plot_name,".pdf"),
          width = 13.33, 
          height = 7.5)
    }
  }else{
    if(!is.na(outname)){
      pdf(paste0(outname, "/", plot_name,".pdf"),
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
  if(!is.na(samplename)){
    if(!is.na(outname)){
      write_csv(data, paste0(outname, "/", samplename, "_", file_name,".csv"))
    }else{
      write_csv(data, paste0(samplename, "_", file_name,".csv"))
    }
  }else{
    if(!is.na(outname)){
      write_csv(data, paste0(outname, "/", file_name,".csv"))
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
if(!is.na(biogenicname) & length(biogenic) > 0){
  message("transforming biogenic data...")
  biogenic_polygon <- image_to_polygon(biogenic, drop_and_fill = F)
  biogenic_raster <- image_to_polygon(biogenic, polygonize = F)
}

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


# element point process models --------------------------------------------
if(start_modeling){
  message("Modeling cell distribution as a function of rock chemistry...")
  
  #generate list of all possible element combinations
  element_combos = list()
  for(element in 1:n_element_combos){
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
    element_models[[i]] <- paste(unlist(element_combos[i]), collapse = "+")
  }
  
  #initialize empty lists for storing model outputs
  element_covariate_models <- list()
  element_model_results <- data.frame(model_ID = as.numeric(), model = as.character(), Npar=as.numeric(), name = as.character(), value = as.numeric()) 
  
  #create the null model for comparison
  ppm0 <- ppm(cell_centroids_ppp ~ 1) 
  
  model_ID = 1
  element_modeler <- function(model){
    message(paste0("Modeling ", model, "..."))
    ppm1 <- ppm(as.formula(paste("cell_centroids_ppp ~", model))) #create a model with covariates
    anova_lrt <- tidy(anova(ppm0,ppm1, test="LRT")) %>% 
      pivot_longer(cols = df:p.value, values_to = "value") %>%
      mutate(model_ID = model_ID, 
             model = model,
             n_elements = str_count(model, "[+]") + 1)
    assign("element_model_results", element_model_results %>%
             bind_rows(anova_lrt), envir = .GlobalEnv)
    assign("model_ID", model_ID + 1, envir = .GlobalEnv)
    element_covariate_models <- append(element_covariate_models, list(ppm1)) #store each model in the element_covariate_models list
  }
  
  element_covariate_models <- lapply(element_models, element_modeler)
  
  
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
  
  saveRDS(best_model_data, file = paste0(outname, "/", samplename,"_best_models.Rds"))
  
  write_data(significant_models, "significant_models")
  
  best_model_elements <- unlist(str_split(best_model$model, "[+]"))
  
  #generate plots 
  pacman::p_load(grid, ggplotlify)
  
  ppm1 <- best_model_data[[1]][[1]]
  
  predict_plot <- as_grob(~plot(predict(ppm1, type = "trend")))
  SE <- predict(ppm1, type = "se")
  SE_plot <- as_grob(~plot(SE, main = "standard error of fitted intensity"))
  
  E_ppm0 <- envelope(ppm0, Linhom, nsim = 159, rank = 1)
  ppm0_plot <- as_grob(~plot(E_ppm0,  main = paste0(samplename, " ppm0")))
  
  E_ppm1 <- envelope(ppm1, Linhom, nsim = 159, rank = 1)
  ppm1_plot <- as_grob(~plot(E_ppm1,  main = paste0(samplename, " ppm1")))
  
  top_row <- gridExtra::grid.arrange(predict_plot,SE_plot, ncol = 2)
  bottom_row <- gridExtra::grid.arrange(ppm0_plot,ppm1_plot, ncol = 2)
  
  arranged_plots <- gridExtra::grid.arrange(top_row, bottom_row, nrow = 2, heights = c(1, 2))
  
  #plot_PDF(grid.draw(arranged_plots), "_ppm1")
  ggsave(paste0(outname, "/", samplename, "_ppm1",".pdf"), 
         grid::grid.draw(arranged_plots), width = 13.33, height = 7.5, units = "in")
  
}

# generate element map with cells  ----------------------------------------


if(plot_maps){
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
  
  if(!is.na(biogenicname) & length(biogenic) > 0){
    density_plot <- density_plot + 
      coord_fixed() +
      ggnewscale::new_scale_color() +
      geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
      scale_color_viridis_c() +
      geom_sf(data = cells_polygon, inherit.aes = F, fill = "#FF00DB", lwd = 0) +
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
      geom_sf(data = cells_polygon, inherit.aes = F, fill = "#FF00DB", lwd = 0) +
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

message("...complete.")
