#load dependencies 
pacman::p_load(shiny, raster, sf, tidyverse)


# import data -------------------------------------------------------------
setwd("/Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO3/D3T13exp_Dec2019_Poorman")
message("importing and aggregating brick...")
xray_brick <- aggregate(brick("D3T13exp_Dec2019_transect_brick.grd"), 50)
SEM <- aggregate(brick("D3T13exp_Dec2019_SEM_pano.tif"), 50)

image_to_polygon <- function(image_path, to_raster = TRUE, pres_abs = TRUE, drop_and_fill = TRUE, polygonize = TRUE, equalizer = TRUE){
  if(equalizer == T){
    message("reading image...")
    image <- image_path %>% image_read() %>% image_quantize(colorspace = 'gray')
    image
    if(to_raster == T){
      message("rasterizing...")
      temp_file <- tempfile()
      image_write(image, path = temp_file, format = 'tiff')
      image_raster_noThresh <- raster(temp_file) 
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
  }else{
    image_read(image_path)
  }
}

cells_polygon <- image_to_polygon("cells.tif", drop_and_fill = F) %>%
  st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")

cell_centroids <- st_centroid(cells_polygon)

cell_centroids_coords <- cell_centroids %>%
  st_coordinates() %>%
  as.data.frame() 


biogenic_polygon <- image_to_polygon("biogenic.tif", drop_and_fill = F) %>%
  st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")


# set color dictionary ----------------------------------------------------

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


# # Define UI for application ---------------------------------------------
ui <- fluidPage(
   
   # Application title
   titlePanel("Element Map"),
   
   # Sidebar with checkboxes for plot features 
   sidebarLayout(
      sidebarPanel(
         checkboxGroupInput("elements", "Choose elements", 
                           choices  = names(xray_brick),
                           selected = c("Fe", "S")),
         checkboxGroupInput("features", "Choose features", 
                            choices  = c("SEM", "cells", "biogenic features", "kernel density"),
                            selected = c("SEM", "kernel density")), 
         downloadButton('down',"download plot"),
         width = 2
      ),
      # Show a plot of the elements
      mainPanel(
         plotOutput("elementPlot"),
         width = 10
         #textOutput("selected_elements")
      )
   )
)


# # Define server logic  --------------------------------------------------
server <- function(input, output) {

# output element plot -----------------------------------------------------
   output$elementPlot <- renderPlot({
     if("SEM" %in% input$features){
       p <-  ggplot() + 
         geom_raster(as.data.frame(SEM, xy=T), mapping=aes(x=x, y=y, fill = layer)) +
         scale_fill_gradient(low = "black",  high = "white") +
         ggnewscale::new_scale_fill()
     }else{
       p <-  ggplot()
     }
     
     
     for(i in input$elements){ 
       p <- p+geom_raster(as.data.frame(xray_brick[[i]], xy=T), mapping=aes_string(x="x", y="y", fill = i, alpha = i)) +
         scale_fill_gradient(low = NA, high = element_colors[i]) + 
         ggnewscale::new_scale_fill()
     }
     
     if("kernel density" %in% input$features){
       p <- p + 
         coord_fixed() +
         ggnewscale::new_scale_color() +
         geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
         scale_color_viridis_c()
     }else{
       p <- p + 
         coord_fixed()
     }
     if("cells" %in% input$features){
       p <- p +
         geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
         coord_sf(datum = NA) 
     }else{
       p <- p 
     }
     if("biogenic features" %in% input$features){
       p <- p +
         geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0)  +
         coord_sf(datum = NA) 
     }else{
       p <- p
     }
     
    p + theme(axis.title = element_blank(),
             axis.text = element_blank(),
             legend.position = "none")
   })
   

# download plot -----------------------------------------------------------
   output$down <- downloadHandler(
     filename = "element_plot.pdf",
     content = function(file) {
       pdf(file, width = 13.33, 
           height = 7.5)
       #xyplot(df[,2]~df[,1],df(),xlim=c(0,10),ylim=c(0,100),type = "b")
       
       # you have to print the plot so that you can open pdf file
         if("SEM" %in% input$features){
           p <-  ggplot() + 
             geom_raster(as.data.frame(SEM, xy=T), mapping=aes(x=x, y=y, fill = layer)) +
             scale_fill_gradient(low = "black",  high = "white") +
             ggnewscale::new_scale_fill()
         }else{
           p <-  ggplot()
         }
         
         
         for(i in input$elements){ 
           p <- p+geom_raster(as.data.frame(xray_brick[[i]], xy=T), mapping=aes_string(x="x", y="y", fill = i, alpha = i)) +
             scale_fill_gradient(low = NA, high = element_colors[i]) + 
             ggnewscale::new_scale_fill()
         }
         
         if("kernel density" %in% input$features){
           p <- p + 
             coord_fixed() +
             ggnewscale::new_scale_color() +
             geom_density_2d(data=cell_centroids_coords, inherit.aes = F, mapping = aes(X, Y, col = stat(level)/max(stat(level)))) +
             scale_color_viridis_c()
         }else{
           p <- p + 
             coord_fixed()
         }
         if("cells" %in% input$features){
           p <- p +
             geom_sf(data = cells_polygon, inherit.aes = F, fill = "#39ff14", lwd = 0) +
             coord_sf(datum = NA) 
         }else{
           p <- p 
         }
         if("biogenic features" %in% input$features){
           p <- p +
             geom_sf(data = biogenic_polygon, inherit.aes = F, fill = "cyan", lwd = 0)  +
             coord_sf(datum = NA) 
         }else{
           p <- p
         }
         
         p <- p + theme(axis.title = element_blank(),
                   axis.text = element_blank(),
                   legend.position = "none")
         print(p)
       dev.off()
     }
   )

   # output$selected_elements <- renderText({
   #   input$elements
   # })
}


# # Run the application  --------------------------------------------------
shinyApp(ui = ui, server = server)

