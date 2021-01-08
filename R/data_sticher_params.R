require(pacman)
pacman::p_load(tidyverse)

file_path <- "/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data"

sample_names <- c("D1T1exp", "D1T2rep", "D1T3rep", "D1T4exp", "D1T5exp", "D1T6rep", "D1T7rep", "D1T8exp",
                  "D3T13exp", "D3T14rep", "D3T15rep", "D3T16exp", "D3T17exp", "D3T18exp", "D3T19rep", "D3T20rep")

stitch_overview <- F

stitch_base_SEM <- F


for(i in 1:length(sample_names)){
  samplename <- sample_names[i]
  site <- if_else(str_detect(samplename, "D1"), "DeMMO1", "DeMMO3")
  directories <- list.dirs(paste(file_path, site, sep = "/"), recursive = F)
  filename <- directories[str_detect(directories, samplename)] #need to detect xray_data
  outname <- paste0("../data/stitched/", sample_names[i])
  source("data_stitcher.R")
}