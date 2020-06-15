pacman::p_load(tidyverse, raster)

directories <- list.dirs("/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data/Dec2019/DeMMO2019_lowMagXEDS", full.names = T , recursive =T)
directories <- directories[2:length(directories)]

files <- list.files(directories, full.names = T)

read_files <- function(file){
  file %>% 
    raster() %>% 
    as.data.frame(xy=T)
}

file_list = lapply(files, read_files)
data <- reduce(file_list, bind_cols) 


test <- str_split(files[1], " ")[[1]][3]

coupon <- str_extract(test, "[^_]+")
scan_id <- str_extract(test, "(?<=_)(.*)(?=[.]tif)")
