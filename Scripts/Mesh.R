library(INLA)
library(dplyr)

setwd("~/Documents/Master_Tesis/INLA_models-master")

#setwd('~/Dropbox/Red lobster project/Data/')
data = read.csv(file = 'datasort_spring_1.csv')


#data = data %>% dplyr::filter(Lat_M < -6.9 & Lat_M > -7)

data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)

#data_1 = data_1 %>% filter(Year>=2014)

data_1 = na.omit(data_1)
data_1$MUN<-ifelse(data_1$MUN>0,1,0)
data_1$BAT<-(data_1$BAT*-1)/1000

#data_1 = data_1 %>% filter(Year>=2004)

#data_1t%>%group_by(Year)%>%summarise(length(which(MUN==0))/length(MUN))
x_loc = cbind(data_1$Lon_M, data_1$Lat_M)
spdat <- data_1
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
# Project data
x_loc_utm <- spTransform(spdat, kmproj)

#Loc = x_loc_utm@coords
#D = dist(Loc)

#hist(D, freq = T, main = '', xlab = 'Distance between points (km)',ylab = 'Frequency')

#plot(x = sort(D), y =(1:length(D))/length(D), type = 'l', xlab = 'Distance between points (km)')


#library(rgdal)
#library(lattice)
#library(ggplot2)
#library(raster)

#setwd("~/Documents/Master_Tesis/mesh ")

#DSN = 'Perfil_peru.shp'
#DSN = 'crop_peru.shp'

#ShapeF.geo = readOGR(dsn = DSN, layer = 'Perfil_peru')

#ShapeF.geo = readOGR(dsn = DSN, layer = 'crop_peru')

#summary(ShapeF.geo)

#proj4string(ShapeF.geo) <- "+init=epsg:4326"
#plot(ShapeF.geo, axes=T)

#sub <- crop(ShapeF.geo, extent(-90, -70, -18, -6.5))
#shapefile(sub, 'crop_peru.shp', overwrite=T)

#ShapeF.utm <- spTransform(ShapeF.geo,kmproj)
#plot(ShapeF.utm, axes=T)

#Peru_df = fortify(ShapeF.utm)
#head(Peru_df)

#coastline = Peru_df[,c('long','lat')]
#plot(x = coastline$long,y=coastline$lat,type='l')

#N = nrow(coastline)
#coast.rev = coastline[N:1, c('long','lat')]
#plot(x = coast.rev$long,y=coast.rev$lat,type='l')

# create basic function
#outerBuffer<-function(x, dist){
#  buff<-buffer(x, width = dist - 1, dissolve = T)
#  e<-erase(buff,x)
#  return(e)
#}

# create data
#p1 = ShapeF.utm

# apply function, create only outside buffer
#a2<-outerBuffer(p1,500)

#saveRDS(a2,'extern_polygon_peru.rds')

# tradaaa ! :)
#plot(a2, col = "black", main= "Outer buffer only")

#meshbuilder()

boundary.loc <- as(x_loc_utm, "SpatialPoints")
boundary <- list(
  inla.nonconvex.hull(coordinates(boundary.loc), 30),
  list(inla.nonconvex.hull(coordinates(boundary.loc), 82)))

## Build the mesh:
mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(30, 62),
                     min.angle=c(30, 25),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=3.6, ## Filter away adjacent points.
                     
                     
# increase cutoff , remove                      
                     
                     offset=c(30, 82)) ## Offset for extra boundaries, if needed.

## Plot the mesh:
plot(mesh, asp=1)


mesh$n  


#setwd('~/Dropbox/Red lobster project/')
setwd("~/Documents/Master_Tesis/INLA_models-master")


saveRDS(mesh,'./Data/mesh_spring_1.rds')

