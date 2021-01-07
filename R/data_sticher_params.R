pacman::p_load(tidyverse)

file_path <- "/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data"


for dir in `find $(pwd) -name 'xray_data'`; 
do sampleID="`basename $(dirname $dir) | cut -f1 -d'_'`";
cd $dir/..;
-f "$dir" -b "$dir/SEM_images" -c "$dir/../coordinates.txt" -z ".tif" -m ".tif" -a "Unknown|SEM|Os" -d "overview" -y "overview" --overview "overview.*tif" -p FALSE -n "$sampleID"; 
done

source("dataStitch.R")