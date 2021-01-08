#load dependencies
pacman::p_load(raster, magick, tidyverse, rasterVis, ggnewscale, Hmisc, cowplot)

#create list all sub-directories within main directory
directories <- list.dirs(filename, full.names = T , recursive =F)

directories <- directories[!str_detect(directories, "Unknown|SEM|Os")]



#set image coordinates
xy <- read_delim(paste0(filename, "../coordinates.txt"), delim = "\t", col_types = cols())
positions <- xy %>%
  select(-x, -y)

message(paste0("Importing data from ", samplename, "..."))

# stitch xray images ------------------------------------------------------

#stitch xray images into panoramas and store in raster brick
xray_brick_list <- list()
xray_data <- list()
for(j in 1:length(directories)){
  path = directories[j]
  files <- list.files(path, full.names = T, pattern = ".tif")
  files <- files[!str_detect(files, "overview")]
  if(length(files) >0){
    xray_data[j] <- str_extract(path, "([^/]+$)")
    message(paste0("Stitching ",xray_data[j], " data (element ", j, " of ", length(directories), ")..."))
    xy_id <- which(positions[[1]] %in% regmatches(files, regexpr("pos.*?\\K-?\\d+", files, perl=TRUE)))
    panorama <- list()
    for(i in 1:length(files)){
      message(paste0("Processing image ", i, " of ", length(files), "..."))
      image_noThresh <- raster(files[i])
      image <- image_noThresh %>% setValues(if_else(values(image_noThresh) > 0, 1, 0))
      image <- aggregate(image, fact=4)
      image_extent <- extent(matrix(c(xy$x[xy_id[i]], xy$x[xy_id[i]] + 1024, xy$y[xy_id[i]], xy$y[xy_id[i]]+704), nrow = 2, ncol = 2, byrow = T))
      image_raster <- setExtent(raster(nrows = 704, ncols = 1024), image_extent, keepres = F)
      values(image_raster) <- values(image)
      panorama[[xy_id[i]]] <- image_raster
    }
      empty_xy_id <- which(!positions[[1]] %in% regmatches(files, regexpr("pos.*?\\K-?\\d+", files, perl=TRUE)))
      if(length(empty_xy_id) > 0){
        for(k in 1:length(empty_xy_id)){
          message(paste0("creating empty raster for missing position: ", empty_xy_id[[k]]))
          empty_raster_extent <- extent(matrix(c(xy$x[empty_xy_id[k]], xy$x[empty_xy_id[k]] + 1024, xy$y[empty_xy_id[k]], xy$y[empty_xy_id[k]]+704), nrow = 2, ncol = 2, byrow = T))
          empty_raster <- setExtent(raster(nrows = 704, ncols = 1024), empty_raster_extent, keepres = F)
          values(empty_raster) <- 0
          panorama[[empty_xy_id[k]]] <- empty_raster
        }
      }
      panorama_merged <- do.call(merge, panorama)
      xray_brick_list[[j]] <- panorama_merged
  }
}

message("Stitching complete. Creating transect x-ray brick...")
xray_brick <- do.call(brick, xray_brick_list)
names(xray_brick) <- xray_data
message("...complete.")


#write the brick
message("Writing transect x-ray brick...")


if(!is.na(outname)){
  dir.create(outname)
  if(!is.na(samplename)){
    out_brick <- writeRaster(xray_brick, paste0(outname, "/", samplename,"_transect_brick.grd"), overwrite=TRUE, format="raster")
    x <- writeRaster(xray_brick, paste0(outname, "/", samplename,"_transect_brick.tif"), overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  }else{
    out_brick <- writeRaster(xray_brick, path = outname, "transect_brick.grd", overwrite=TRUE, format="raster")
    x <- writeRaster(xray_brick, path = outname, "transect_brick.tif", overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  }
}else{
  if(!is.na(samplename)){
    out_brick <- writeRaster(xray_brick, paste0(samplename,"_transect_brick.grd"), overwrite=TRUE, format="raster")
    x <- writeRaster(xray_brick, paste0(samplename,"_transect_brick.tif"), overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  }else{
    out_brick <- writeRaster(xray_brick, "transect_brick.grd", overwrite=TRUE, format="raster")
    x <- writeRaster(xray_brick, "transect_brick.tif", overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  }
}


message("...complete.")

if(stitch_overview){
  message("Stacking overview x-ray brick...")
  overview_brick_list <- list()
  overview_data <- list()

  for(i in 1:length(directories)){
    path = directories[i]
    files <- list.files(path, full.names = T, pattern = "overview.*tif")
    if(length(files) >0){
      overview_data[[i]] <- str_extract(path, "([^/]+$)")
      message(paste0("Stacking ",overview_data[[i]], " data (element ", i, " of ", length(directories), ")..."))
      image_noThresh <- raster(files)
      image <- image_noThresh %>% setValues(if_else(values(image_noThresh) > 0, 1, 0))
      image <- aggregate(image, fact=4)
      image_extent <- extent(matrix(c(0, 1024, 0, 704), nrow = 2, ncol = 2, byrow = T))
      image_raster <- setExtent(raster(nrows = 704, ncols = 1024), image_extent, keepres = F)
      values(image_raster) <- values(image)
      overview_brick_list[[i]] <- image_raster
    }
  }

  message("Stacking complete. Creating overview x-ray brick...")
  overview_brick <- do.call(brick, na.omit(overview_brick_list))
  names(overview_brick) <- unlist(overview_data)

  message("...complete. Writing overview x-ray brick...")
  if(!is.na(outname)){
    dir.create(outname)
    if(!is.na(samplename)){
      out_brick <- writeRaster(overview_brick, paste0(outname, "/", samplename,"_overview_brick.grd"), overwrite=TRUE, format="raster")
      x <- writeRaster(overview_brick, paste0(outname, "/", samplename,"_overview_brick.tif"), overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    }else{
      out_brick <- writeRaster(overview_brick, path = outname, "overview_brick.grd", overwrite=TRUE, format="raster")
      x <- writeRaster(overview_brick, path = outname, "overview_brick.tif", overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    }
  }else{
    if(!is.na(samplename)){
      out_brick <- writeRaster(overview_brick, paste0(samplename,"_overview_brick.grd"), overwrite=TRUE, format="raster")
      x <- writeRaster(overview_brick, paste0(samplename,"_overview_brick.tif"), overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    }else{
      out_brick <- writeRaster(overview_brick, "overview_brick.grd", overwrite=TRUE, format="raster")
      x <- writeRaster(overview_brick, "overview_brick.tif", overwrite=TRUE, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    }
  }
  message("...complete.")
}


#flush everything we don't need from memory
remove(list = c("x", "xray_brick_list", "xray_data", "overview_brick_list", "overview_data", "empty_raster", "empty_raster_extent",
                "i", "j", "k", "path", "positions", "xy_id", 
                "panorama_merged", "panorama", "empty_xy_id", "files", "directories"))


# create base SEM image ---------------------------------------------------
  

SEM_images <- list.files(paste0(filename, "/SEM_images"), full.names = T, pattern = ".tif")


SEM_images  <- SEM_images[!str_detect(SEM_images , "overview")]

SEM_panorama <- list()

message("Stitching SEM images into panorama...")
for(i in 1:length(SEM_images)){
  image <- raster(SEM_images[i])
  image_extent <- extent(matrix(c(xy$x[i], xy$x[i] + 1024, xy$y[i], xy$y[i]+704), nrow = 2, ncol = 2, byrow = T))
  image_raster <- setExtent(raster(nrows = 704, ncols = 1024), image_extent, keepres = F)
  values(image_raster) <- values(image)
  SEM_panorama[[i]] <- image_raster
}
SEM_panorama_merged <- do.call(merge, SEM_panorama)
message("...complete.")

#write the brick 
message("Writing SEM panoramic raster...")


if(!is.na(outname)){
  if(!is.na(samplename)){
    writeRaster(SEM_panorama_merged, paste0(outname, "/", samplename,"_SEM_pano.tif"), overwrite=TRUE, format = "GTiff")
  }else{
    writeRaster(SEM_panorama_merged, path = outname, "SEM_pano.tif", overwrite=TRUE, format = "GTiff")
  }
}else{
  if(!is.na(samplename)){
    writeRaster(SEM_panorama_merged, paste0(samplename, "_SEM_pano.tif"), overwrite=TRUE, format = "GTiff")
  }else{
    writeRaster(SEM_panorama_merged, "SEM_pano.tif", overwrite=TRUE, format = "GTiff")
  }
}


message("...complete.")


# create SEM overview raster ----------------------------------------------

if(stitch_overview)){
  SEM_images <- list.files(paste0(filename, "/SEM_images"), full.names = T, pattern = ".tif")

  message("Creating overview SEM raster image...")
  SEM_images <- SEM_images[str_detect(SEM_images , "overview.*tif")]
  image <- raster(SEM_images)
  image_extent <- extent(matrix(c(0, 1024, 0, 704), nrow = 2, ncol = 2, byrow = T))
  image_raster <- setExtent(raster(nrows = 704, ncols = 1024), image_extent, keepres = F)
  values(image_raster) <- values(image)
  SEM_image <- image_raster
  
  message("Writing overview SEM raster...")
  
  
  if(!is.na(outname)){
    if(!is.na(samplename)){
      writeRaster(SEM_image, paste0(outname, "/", samplename,"_SEM_overview.tif"), overwrite=TRUE, format = "GTiff")
    }else{
      writeRaster(SEM_image, path = outname, "SEM_overview.tif", overwrite=TRUE, format = "GTiff")
    }
  }else{
    if(!is.na(samplename)){
      writeRaster(SEM_image, paste0(samplename, "_SEM_overview.tif"), overwrite=TRUE, format = "GTiff")
    }else{
      writeRaster(SEM_image, "SEM_overviewo.tif", overwrite=TRUE, format = "GTiff")
    }
  }
}
