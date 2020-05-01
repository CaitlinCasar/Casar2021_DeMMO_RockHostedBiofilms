## to run a single job with dataStitchR

Rscript /Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/R/scipts/dataStitchR.R -f "xray_data" -b "xray_data/SEM_images" -c "coordinates.txt" -z ".tif" -m ".tif" -a "Unknown|Os|SEM" -d "overview" -y "overview" --overview "overview.*tif" -p FALSE -n "D1T4exp"

## to run a single job with dataClustR

Rscript /Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/R/scipts/dataClustR.R -f "/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data/DeMMO3/D3T20rep_Dec2019_Yates" --transect "transect.*grd" --overview "overview.*grd" --out "reports" -n "D3T20rep" --base "SEM_pano.tif" -u 4.2 --cores 3 --model_vars 5 --cells "cells.tif" --biogenic "biogenic.tif" -p TRUE



## to run batch job for dataStitchR

cd /Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data

for dir in `find $(pwd) -name 'xray_data'`; 
do sampleID="`basename $(dirname $dir) | cut -f1 -d'_'`";
cd $dir/..;
Rscript /Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/R/scipts/dataStitchR.R -f "$dir" -b "$dir/SEM_images" -c "$dir/../coordinates.txt" -z ".tif" -m ".tif" -a "Unknown|SEM|Os" -d "overview" -y "overview" --overview "overview.*tif" -p FALSE -n "$sampleID"; 
done


## to run batch job for dataClustR

cd /Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data

for dir in `find $(pwd) -name '*Dec2019*'`; 
do sampleID="`basename $dir | cut -f1 -d'_'`";
cd $dir;
Rscript /Users/Caitlin/Desktop/DeMMO_Pubs/DeMMO_NativeRock/DeMMO_NativeRock/R/scipts/dataClustR.R -f "$dir" --transect "transect.*grd" --overview "overview.*grd" --out "reports" -n "$sampleID" --base "SEM_pano.tif" -u 4.2 --cores 3 --model_vars 5 --cells "cells.tif" --biogenic "biogenic.tif" -p TRUE
done



