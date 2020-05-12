library(raster)
library(rgeos)

# example data
x <- raster("../test_images/Fe.tif")
#p <- rasterToPolygons(x, dissolve=TRUE)

### To get the rectangular extent
e <- extent(x)
# coerce to a SpatialPolygons object
p <- as(e, 'SpatialPolygons')  

# calculate the centroid 
c1 <- gCentroid(p)
plot(p)
plot(c1, col='blue', add=TRUE)


# make all values the same. Either do
r <- x > -Inf
# or alternatively
# r <- reclassify(x, cbind(-Inf, Inf, 1))

# convert to polygons (you need to have package 'rgeos' installed for this to work)
pp <- rasterToPolygons(r, dissolve=TRUE, fun=function(x){x<6})

# look at the results
plot(x)
plot(p, lwd=5, border='red', add=TRUE)
plot(pp, lwd=3, border='blue', add=TRUE)




# x -----------------------------------------------------------------------



library(raster)
library(stars)
library(sf)
library(magrittr)

f <- system.file("../test_images/Fe.tif", package="raster")
r <- raster(f)
r[r[] < 750] <- 0
r[r[] >= 750] <- 1

z <- st_as_stars(x) %>% 
  st_as_sf() %>% # this is the raster to polygons part
  st_cast("MULTILINESTRING") # cast the polygons to polylines

plot(z)





# x -----------------------------------------------------------------------
library(xROI)
test_array <- getCLArray("../test_images/Fe_colored.jpg")

f <- system.file(package = 'xROI', "../test_images/Fe_colored.jpg")
cli <- getCL(f)

#generate a polygon 
polygon_coords <- if(interactive()){
  drawPolygon()
}



# x -----------------------------------------------------------------------

library(raster)
library(rgeos)

# get Brasil borders

shp <- getData(country = 'BRA',level=0)

#create binary raster


r <- raster(extent(shp),resolution=c(0.5,0.5))
r[] <- NA # values have to be NA for the buffering

# take centroid of Brasil as center of species presence
cent <- gCentroid(shp)

# set to 1 
r[cellFromXY(r,cent)] <- 1

# buffer presence
r <- buffer(r,width=1000000)

# set rest 0
r[is.na(r)] <- 0

# mask by borders
r <- mask(r,shp)

pol <- rasterToPolygons(r,function(x) x == 1,dissolve=T)

plot(shp)
plot(pol,col='red',add=T)

# adjust tol for smoothness
pol_sm <- gSimplify(pol,tol=0.5)

plot(pol)
lines(pol_sm,col='red',lwd=2)




