library (INLA)
INLA:::inla.dynload.workaround()
library (dplyr)
library(rgeos)
library(sf) # for vector data 
library(raster) # for working with rasters
library(maps) # additional helpful mapping packages
library(tidyr)
source('./Scripts/spde-book-functions.R')

# Input data #
#data_c <-read.csv(file = 'Data/data_complete.csv')
# only distribution of RSL
#data = data_c %>% dplyr::filter(Lat_M < -6.5)
# Only summer (more overlapping)
#data = data %>% filter(season=='Summer')
# Only night (more overlapping)
#data = data %>% filter(day_night=='Night')
##only 6 years
#data =data %>% filter(Year>2012)

data = read.csv(file = 'data_filtered_complete.csv')
data<-data%>%transmute(Lon_M,Lat_M,MUN,ANC,SST,SALI,GSST, Year,DC,BAT,CHL)
data = na.omit(data)# si no filtras antes las variables que quieres estas eliminado el 60% de tus datos

x_loc = cbind(data$Lon_M, data$Lat_M)
spdat <- data
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
# Project data
x_loc_utm <- spTransform(spdat, kmproj)

# mesh building 
data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)
data_1<-data_1%>%filter(Year>2010)
#data_1 = na.omit(data_1)
data_1$MUN<-ifelse(data_1$MUN>0,1,0)
data_1$BAT<-(data_1$BAT*-1)/1000
RNGkind(sample.kind = 'Rounding')
set.seed(1234)
sample_size = floor(0.2*nrow(data_1))
samplet<- sample(seq_len(nrow(data_1)),size = sample_size)
data_1t<-data_1[samplet,]
data_1p<-data_1[-samplet,]

#write.csv(data_1t,'~/Dropbox/Red lobster project/Data/data_1t.csv')
#data_1t = read.csv('~/Dropbox/Red lobster project/Data/data_1t.csv')
data_1t<-data_1t[order(data_1t$Year),]
table(data_1t$Year)
x_loc = cbind(data_1t$Lon_M, data_1t$Lat_M)
spdat <- data_1t
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
# Project data
x_loc_utm <- spTransform(spdat, kmproj)

plot(x_loc_utm)

#meshbuilder()

# ## Build the mesh:
# boundary.loc <- as(x_loc_utm, "SpatialPoints")
# boundary <- list(
#   inla.nonconvex.hull(coordinates(boundary.loc), 40),
#   inla.nonconvex.hull(coordinates(boundary.loc), 78))
# 
# ## Build the mesh:
# mesh <- inla.mesh.2d(boundary=boundary,
#                      max.edge=c(20, 62),
#                      min.angle=c(30, 25),
#                      cutoff=1.8, ## Filter away adjacent points.
#                      offset=c(20, 40)) ## Offset for extra boundaries, if needed.
# 
#saveRDS(mesh,'mesh.rds')

mesh = readRDS('mesh.rds') 

## Plot the mesh:
plot(mesh,asp=1)

mesh$n

#2.1 Build the spde
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(60, 0.1), # is (range0,Prange) P(range < range0) = Prange ## 25 km instead of 40 
                            prior.sigma = c(1, 0.01)) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
# 
#Prior r0 better not to favour ranges that are smaller than the resolution of the mesh
#Prior sigma how variable is the field from a point to the next

#2.2 Define the time effects,NEW
t<-8#18
mesh.t1 = inla.mesh.1d(1:8)
data_1t$Year<-as.numeric(factor(data_1t$Year))
#3. Projector matrix
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm, group.mesh = mesh.t1,#anadi group.mesh
                         group = data_1t$Year)#NEW

dim(Aspte)

#4.Stack data
#Xm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_1t)
N = nrow(data_1t)
#X = data.frame(Intercept = rep(1,N), DC=Xm[,1], SST=Xm[,2], SALI=Xm[,3], GSST=Xm[,4], BAT=Xm[,5],CHL=Xm[,6]) 

index.st <- inla.spde.make.index('i', n.spde = mesh$n, n.group = 8)

#NEW INDEX SPATIO TEMPORAL
stk.dat.st <- inla.stack(tag='dat', data = list(z = data_1t$MUN), 
                         A = list(Aspte,1), effects = list(index.st,list(Intercept=rep(1,N),
                         DC=data_1t$DC,BAT=data_1t$BAT,SST=data_1t$SST,SALI=data_1t$SALI,
                         CHL=data_1t$CHL,GSST=data_1t$GSST))) #

#NEW INDEX SPATIO TEMPORAL  stk.dat.st <- inla.stack(data = list(z = data06_13t$Anch),    A = list(Aspte,1,1),   effects = list(index.st, #NEW   Intercept = 1,X = X),tag = 'dat') 


# #4.1 stack for prediction
# #A projector for prediction
time_mesh<-sort(rep(seq(1:10),mesh$n))

mesh.loc =cbind(rep(mesh$loc[,1],10),rep(mesh$loc[,2],10),sort(rep(seq(1:10),mesh$n)))
#shoreline<-read.csv("Data/shoreline1.csv",sep=",",header = TRUE)
#coordinates(shoreline) <- ~ X + Y
#proj4string(shoreline) <- "+init=epsg:4326"
#shoreline_utm = spTransform(shoreline, kmproj)
meshcoord <-data.frame(mesh$loc[,1],mesh$loc[,2])
names(meshcoord)<-c('utmx','utmy')

#write.csv(meshcoord,'meshcoordutm_100.csv')

coordinates(meshcoord) <- ~ utmx + utmy
proj4string(meshcoord) <- "+proj=utm +zone=18S ellps=WGS84 +units=km +ellps=WGS84"

#DCmesh <- rep(apply(spDists(shoreline_utm,meshcoord,longlat = FALSE), 2, min),18)#only run the first time

Aprd <- inla.spde.make.A(mesh,loc=mesh.loc,n.group=mesh.t1$n,group=time_mesh)#

dim(Aprd)

data_pred = read.csv('var_pred_c.csv')
data_pred1<-data_pred[64176:67950,]
data_at<-data.frame(x=meshcoord@coords[,1],y=meshcoord@coords[,2],data_pred1)
ggplot(aes(x,y),data=data_at)+geom_point(aes(colour=CHL))
#4.Stack data
#Xpm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_pred)
N_pred = nrow(data_pred1)
#Xp = data.frame(Intercept = rep(1,N_pred), DC=Xpm[,1], SST=Xpm[,2], SALI=Xpm[,3], GSST=Xpm[,4], BAT=Xpm[,5],CHL=Xpm[,6]) 


# # #Mesh location
# #stack prediction
stk.prd.st <- inla.stack(data = list(z = NA),
                         A = list(Aprd,1), 
                         effects = list(index.st,#,mesh.loc
                         list(Intercept = rep(1,N_pred),
                         DC=data_pred1$DC,SST=data_pred1$SST,SALI=data_pred1$SALI,CHL=data_pred1$CHL)),tag = 'prd') #

#4.2 stack all
#stk.all.st <- inla.stack(stk.dat.st)
stk.all.st <- inla.stack(stk.dat.st, stk.prd.st)#join all in stack

#5.Formula
f.DCst <- z ~ -1 + Intercept + DC + BAT+SST+ SALI +CHL+GSST+
  f(i, model = spde,   group = i.group, control.group = list(model = 'ar1',
    hyper=list(theta=list(prior='pc.cor1', param=c(0.3, 0.7)))))#NEW

#parameters of the pc prior for correlation are u=sqrt(3)/2 and alpha = 3/4 
#P(r>r0)=p,Prob(cor>u)=alpha
#6.Fitting the model
r.DCst <- inla(f.DCst, family = 'binomial', 
               control.family=list(link='logit'),
               data = inla.stack.data(stk.dat.st), 
               control.predictor = list(A = inla.stack.A(stk.dat.st), link = 1,#the link 1 is to indicate that fitted values are required
                                        compute=TRUE),
               control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),verbose=TRUE)
summary(r.DCst)
saveRDS(r.DCst,'bin_1119.rds')
#Plot the posterior
par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.DCst$marginals.fix[[1]], type = 'l', #only 75 values in marginal fixed always
     xlab = 'Intercept', ylab = 'Density') 
plot(r.DCst$marginals.fix[[2]], type = 'l', 
     xlab = 'Ocean distance coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[3]], type = 'l', 
     xlab = 'Batimetry  coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[4]], type = 'l', 
     xlab = 'SST  coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[5]], type = 'l', 
     xlab = 'Salinity coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[6]], type = 'l', 
     xlab = 'Chl', ylab = 'Density') 
plot(r.DCst$marginals.fix[[7]], type = 'l', 
     xlab = 'GSST', ylab = 'Density') 
par(mfrow = c(3, 1), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
for (j in 1:3) 
  plot(r.DCst$marginals.hy[[j]], type = 'l', 
       ylab = 'Density', xlab = names(r.DCst$marginals.hyperpar)[j])

r.DCst$summary.random$i$mean
plot(r.DCst$summary.random$i$mean)#Plot of the field,represent the spatial structure
#7.Prediction on a grid
## ----GRID for projection of field --------------------------------------------------
#GRID 10km###
stepsize <- 10 #10km
poly.only = readRDS('Data/poly.only')
poly.grid<-spatialEco::remove.holes(poly.only)
poly.gridutm = spTransform(poly.grid, kmproj)# Project data
x.range <- diff(poly.gridutm@bbox[1,])
y.range <- diff(poly.gridutm@bbox[1,])
nxy <- round(c(x.range, y.range) / stepsize)
projgrid <- inla.mesh.projector(mesh, xlim = poly.gridutm@bbox[1,], 
                                ylim = poly.gridutm@bbox[2,], dims = nxy)
library(splancs)
xy.in <- inout(projgrid$lattice$loc, as.matrix(poly.gridutm@polygons[[1]]@Polygons[[1]]@coords))#Test points for inclusion in a polygon
table(xy.in)
mu.st1 <- lapply(1:mesh$n, function(j) {
  idx <- 1:spde$n.spde + (j - 1) * spde$n.spde
  r <- inla.mesh.project(projgrid,field = r.DCst$summary.ran$i$mean[idx])
  r[!xy.in] <- NA
  return(r)
})
years<-c('2011','2012','2013','2015','2016','2017','2018','2019')
par(mfrow = c(4, 2), mar = c(0, 0, 1, 0))
zlm1 <- range(unlist(mu.st1), na.rm = TRUE)
# identify wich time location is near each knot
for (j in 1:8) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st1[[j]]), 
                  zlim = zlm1, main = years[j])
}

#Prediction of the response by computation of posterior distribution
## ------------------------------------------------------------------------
source('./Scripts/spde-book-functions.R')
# #Matrix for prediction
index_predic <- inla.stack.index(stk.all.st,"prd")$data
stpred <- matrix(r.DCst$summary.fitted.values$mean[index_predic], 
                 spde$n.spde)

par(mfrow = c(6, 3), mar = c(1, 1, 1, 1))
for(i in 1:6){
  meanproj <- inla.mesh.project(projgrid,field=stpred[, i])
  meanproj[!xy.in] <- NA
  print(mean(meanproj,na.rm=TRUE))
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = meanproj))
}

###New part

summary(proof3)
#Plot the posterior
par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(proof3$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept', ylab = 'Density')
plot(proof3$marginals.fix[[2]], type = 'l', 
     xlab = 'DC', ylab = 'Density')
plot(proof3$marginals.fix[[3]], type = 'l', 
     xlab = 'SST', ylab = 'Density') 
plot(proof3$marginals.fix[[4]], type = 'l', 
     xlab = 'SAL', ylab = 'Density') 
plot(proof3$marginals.fix[[5]], type = 'l', 
     xlab = 'Bat', ylab = 'Density') 
plot(proof3$marginals.fix[[6]], type = 'l', 
     xlab = 'GSST', ylab = 'Density') 
plot(proof3$marginals.fix[[7]], type = 'l', 
     xlab = 'Chla', ylab = 'Density') 
par(mfrow = c(3, 3), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
for (j in 1:3) 
  plot(proof3$marginals.hy[[j]], type = 'l', xlab = names(proof3$marginals.hyperpar)[j])

