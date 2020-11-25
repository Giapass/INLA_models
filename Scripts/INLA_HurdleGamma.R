# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
###########################################################
# Load packages
library (INLA)
library (dplyr)
library(tidyr)
library(rgeos)
library(sf) # for vector data 
library(raster) # for working with rasters
library(maps) # additional helpful mapping packages
source('F:/Postdoc Hamburg/Master Pleuroncodes/spde-book-files/R/spde-book-functions.R')
###Input data#########
data06_13 <-readRDS(file = "./Data/data06_13.rds")
data06_13<-data06_13%>%dplyr::filter(Lon>-84.47361)%>%drop_na()
data06_13$z <- (data06_13$Anch>0) + 0
data06_13$y<- ifelse(data06_13$Anch>0, data06_13$Anch, NA)
set.seed(123)
sample_size = floor(0.2*nrow(data06_13))
samplet<- sample(seq_len(nrow(data06_13)),size = sample_size)
data06_13t<-data06_13[samplet,]
data06_13t<-data06_13t[order(data06_13t$Year),]
data06_13p<-data06_13[-samplet,]
x_loc = cbind(data06_13t$Lon, data06_13t$Lat)
spdat <- data06_13t
coordinates(spdat) <- ~ Lon + Lat
proj4string(spdat) <- "+init=epsg:4326"
#####Priors######################
prange <- c(40, 0.01)#spde
psigma<-c(1, 0.3)#spde
rhoprior <- list(theta = list(prior = 'pccor1',param = c(0.2, 0.7)))#rho
bprior <- list(prior = 'gaussian', param = c(0,1))#betaprior, shared
pcgprior <- list(prior = 'pc.gamma', param = 1)#gamma prior
##################################
#1.Define the boundaries
library(raster)
filename <- "./Data/peru_extend.shp"
peru_extend<- shapefile(filename)
# #Plot Polygons
pl.grid <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(-84.47361,-84.47361,-70.37,-70.37), 
        c(-3.3921,-18.41958,-18.41958,-3.3921)),
  FALSE)), '0')), proj4string = CRS(proj4string(peru_extend)))
poly.only<- gDifference(pl.grid, peru_extend)
# writeOGR(poly.only, dsn = '.', layer = 'poly', driver = "ESRI Shapefile")
# # #Plot data
plot(poly.only);points(data06_13t$Lon, data06_13t$Lat,cex=0.1,col=2)
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
# Project data
poly.onlyutm = spTransform(poly.only, kmproj)
pl.sel = spTransform(pl.grid, kmproj)
## ---------------------Define the mesh---------------------------------------------------
mesh <- inla.mesh.2d(boundary = poly.onlyutm,
                     max.edge = c(3,9) * 10,
                     cutoff = 2,# Filter away adjacent points
                     offset = c(30, 150))#
plot(mesh);mesh$n#number of vertices
dim(mesh$graph$vv)#neighborhood
mesh$idx$loc#no one because the mesh is not build according to the points
x_loc_utm <- spTransform(spdat, kmproj)
#2.1 Build the spde
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = prange, # is (range0,Prange) P(range < range0) = Prange
                            prior.sigma = psigma) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
#Prior r0 better not to favour ranges that are smaller than the resolution of the mesh
#Prior sigma how variable is the field from a point to the next
#2.2 Define the time effects,NEW
t<-6
mesh.t1 = inla.mesh.1d(1:6)
data06_13t$time<-as.numeric(factor(data06_13t$Year))
#3. Projector matrix#######
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm,
                         group.mesh = mesh.t1, group = data06_13t$time)#NEW
dim(Aspte)==c(length(data06_13t$Year),(6*spde$n.spde))#to confirm the correct dim of matrix A
#4.Stack data##############
#4.1 Make index
index.stz <- inla.spde.make.index('i_z', n.spde = spde$n.spde,
                                 n.group = mesh.t1$n)#NEW INDEX SPATIO TEMPORAL
index.stc <- inla.spde.make.index('i_c', n.spde = spde$n.spde,
                                  n.group = mesh.t1$n)#NEW INDEX SPATIO TEMPORAL
index.sty <- inla.spde.make.index('i_y', n.spde = spde$n.spde,
                                  n.group = mesh.t1$n)#NEW INDEX SPATIO TEMPORAL
#4.2 Stack estimation
stk.z.est <- inla.stack(data = list(Y = cbind(data06_13t$z,NA),link=1),#NEW 
                         A = list(Aspte,1),
                         effects = list(index.stz, 
                                        data.frame(z.Intercept =rep(1,length(data06_13t$DC)),oceanDist_z = data06_13t$DC)),tag = 'zdat') 
stk.y.est <- inla.stack(data = list(Y = cbind(NA,data06_13t$y),link=2),#NEW 
                        A = list(Aspte,1),
                        effects = list(c(index.stc,index.sty), 
                                       data.frame(y.Intercept = rep(1,length(data06_13t$DC)),oceanDist_y = data06_13t$DC)),tag = 'ydat') 
#4.3 Stack prediction
time_mesh<-sort(rep(seq(1:6),mesh$n))
mesh.loc =cbind(rep(mesh$loc[,1],6),rep(mesh$loc[,2],6),sort(rep(seq(1:6),mesh$n)))
shoreline<-read.csv("F:/Postdoc Hamburg/INLA/shoreline1.csv",sep=",",header = TRUE)
coordinates(shoreline) <- ~ X + Y
proj4string(shoreline) <- "+init=epsg:4326"
shoreline_utm = spTransform(shoreline, kmproj)
meshcoord <-data.frame(mesh$loc[,1],mesh$loc[,2])
names(meshcoord)<-c('utmx','utmy')
coordinates(meshcoord) <- ~ utmx + utmy
proj4string(meshcoord) <- "+proj=utm +zone=18S ellps=WGS84 +units=km +ellps=WGS84"
DCmesh <- rep(apply(spDists(shoreline_utm,meshcoord,longlat = FALSE), 2, min),6)#only run the first time
Aprd <- inla.spde.make.A(mesh,loc=mesh.loc,n.group=mesh.t1$n,group=time_mesh)#
stk.z.prd <- inla.stack(data = list(Y = matrix(NA, ncol(Aprd), 2), link = 1),
                         A = list(Aprd,1), 
                         effects = list(index.stz,
                                        data.frame(z.Intercept=rep(1,length(DCmesh)),oceanDist_z=DCmesh)),tag = 'zprd') 
stk.y.prd <- inla.stack(data = list(Y = matrix(NA, ncol(Aprd), 2), link = 2),
                         A = list(Aprd,1), 
                         effects = list(c(index.stc,index.sty),
                                        data.frame(y.Intercept=rep(1,length(DCmesh)),oceanDist_y=DCmesh)),tag = 'yprd') 
#4.4 Stack all
stk.all.hurdle <- inla.stack(stk.z.est,stk.y.est, stk.z.prd, stk.y.prd)

#5.Formula
f.DCjoint <- Y ~ -1 + z.Intercept +y.Intercept+oceanDist_z+oceanDist_y+
  f(i_z,model=spde, group=i_z.group, control.group = list(model = 'ar1', hyper = rhoprior))+
  f(i_c, copy = "i_z", group = i_c.group, hyper = list(theta = bprior)) + 
  f(i_y, model = spde, group = i_y.group, control.group =list(model = 'ar1', hyper = rhoprior))  

#6.Fitting the model
r.DCjoint <- inla(f.DCjoint, family =c("binomial", "gamma"), 
               data = inla.stack.data(stk.all.hurdle), 
               control.family =list(list(), list(hyper = list(theta = pcgprior))), 
               control.predictor = list(A = inla.stack.A(stk.all.hurdle), link = 1,compute=TRUE),
               control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),verbose=TRUE)
summary(r.DCjoint)
#Plot the posterior
par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.DCjoint$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept.Z', ylab = 'Density')
plot(r.DCjoint$marginals.fix[[2]], type = 'l', 
     xlab = 'Intercept.y', ylab = 'Density')
plot(r.DCjoint$marginals.fix[[3]], type = 'l', 
     xlab = 'Ocean distance coefficient', ylab = 'Density') 
for (j in 1:4) 
  plot(r.DCjoint$marginals.hy[[j]], type = 'l', 
       ylab = 'Density', xlab = names(r.DCjoint$marginals.hyperpar)[j])
#7.Prediction on a grid
## ----GRID for projection of field --------------------------------------------------
#GRID 10km###
stepsize <- 10 #10km
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
  r <- inla.mesh.project(projgrid,field = r.DCjoint$summary.ran$i$mean[idx])
  r[!xy.in] <- NA
  return(r)
})
par(mfrow = c(4, 4), mar = c(0, 0, 1, 0))
zlm1 <- range(unlist(mu.st1), na.rm = TRUE)
# identify wich time location is near each knot
for (j in 1:6) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st1[[j]]), 
                  zlim = zlm1, main = paste0("Mesh timepoint: ", j))
}

#Prediction of the response by computation of posterior distribution
## ------------------------------------------------------------------------
source('F:/Postdoc Hamburg/Master Pleuroncodes/spde-book-files/R/spde-book-functions.R')
# #Matrix for prediction
idx.zdat <- inla.stack.index(stk.all.hurdle, 'zdat')$data
idx.ydat <- inla.stack.index(stk.all.hurdle, 'ydat')$data
idx.zprd <- inla.stack.index(stk.all.hurdle, 'zprd')$data
idx.yprd <- inla.stack.index(stk.all.hurdle, 'yprd')$data


stprdz <- matrix(r.DCjoint$summary.fitted.values$mean[idx.zprd], spde$n.spde)
stprdy <- matrix(r.DCjoint$summary.fitted.values$mean[idx.yprd], spde$n.spde)

par(mfrow = c(3, 2), mar = c(0, 0, 0, 0))
for(i in 1:6){
  meanproj <- inla.mesh.project(projgrid,field=stprdz[, i])
  meanproj[!xy.in] <- NA
  print(mean(meanproj,na.rm=TRUE))
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = meanproj))
}

par(mfrow = c(3, 2), mar = c(0, 0, 0, 0))
for(i in 1:6){
  meanproj <- inla.mesh.project(projgrid,field=stprdy[, i])
  meanproj[!xy.in] <- NA
  print(mean(meanproj,na.rm=TRUE))
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = meanproj))
}




