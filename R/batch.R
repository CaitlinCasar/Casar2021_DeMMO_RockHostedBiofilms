require(pacman)
pacman::p_load(tidyverse)


# set file path and samples to be processed -------------------------------

sample_names <- c("D1T1exp", "D1T2rep", "D1T3rep", "D1T4exp", "D1T5exp", "D1T6rep", "D1T7rep", "D1T8exp",
                  "D3T13exp", "D3T14rep", "D3T15rep", "D3T16exp", "D3T17exp", "D3T18exp", "D3T19rep", "D3T20rep")


file_path <- "/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data"


# stitch image and X-ray data ---------------------------------------------

stitch_overview <- F


for(i in 1:length(sample_names)){
  samplename <- sample_names[i]
  site <- if_else(str_detect(samplename, "D1"), "DeMMO1", "DeMMO3")
  directories <- list.dirs(paste(file_path, site, sep = "/"), recursive = F)
  filename <- directories[str_detect(directories, samplename)] #need to detect xray_data
  outname <- paste0("../data/stitched/", sample_names[i])
  source("data_stitcher.R")
}


# model data + estimate 2D kernel densities -------------------------------

transectname <- "transect.*grd" 

overviewname <- "overview.*grd"

overview_base <- "SEM_overview.tif"

base_name <- "SEM_pano.tif"

cellsname <- "cells.tif"

biogenicname <- "biogenic.tif"

scale_conversion <- 4.2

start_modeling <- F

plot_maps <- T

n_element_combos <- 5


for(i in 1:length(sample_names)){
  samplename <- sample_names[i]
  site <- if_else(str_detect(samplename, "D1"), "DeMMO1", "DeMMO3")
  directories <- list.dirs(paste(file_path, site, sep = "/"), recursive = F)
  filename <- directories[str_detect(directories, samplename)]
  outname <- paste0("../data/reports/", sample_names[i])
  source("element_modeler.R")
}


# calculate bulk rock compositions ----------------------------------------

source("bulkElementComposition.R")


# asv diversity -----------------------------------------------------------

source("asvs.R")


# summarize data ----------------------------------------------------------

source("summarize_data.R")

