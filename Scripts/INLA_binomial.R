# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
###########################################################
# Load packages
  library (INLA)
  library (dplyr)
  library(rgeos)
  library(sf) # for vector data 
  library(raster) # for working with rasters
  library(maps) # additional helpful mapping packages
library(tidyr)
source('F:/Postdoc Hamburg/Master Pleuroncodes/INLA_models/Scripts/spde-book-functions.R')
###Input data########
data<-readRDS(file = "./Data/data06_13.rds")
data06<-data%>%filter(Year==2006)%>%dplyr::select(Year,Lon,Lat,Anch,DC)%>%drop_na()
data06$Anch<-ifelse(data06$Anch>0,1,0)
set.seed(123)
sample_size = floor(0.7*nrow(data06))
samplet<- sample(seq_len(nrow(data06)),size = sample_size)
data06t<-data06[samplet,]
data06p<-data06[-samplet,]
x_loc = cbind(data06t$Lon, data06t$Lat)
spdat <- data06t
coordinates(spdat) <- ~ Lon + Lat
proj4string(spdat) <- "+init=epsg:4326"
#---Simple binomial model------
model.simple<-inla(Anch~DC,data=data06t,family = "binomial",
control.compute = list(waic=T,dic=T,cpo=T),
control.predictor = list(compute = TRUE, link=1))
summary(model.simple)
plot(model.simple$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept', ylab = 'Density') 
plot(model.simple$marginals.fix[[2]], type = 'l', 
     xlab = 'Lat', ylab = 'Density') 
plot(data06t$DC,model.simple$summary.fitted.values$mean)
#-----Binomial model with spatial effect------------------
#1.The mesh for predictions
  library(raster)
  pe <- getData('GADM', country='PE', level=0)
  br <- getData('GADM', country='BR', level=1)
  br_r <- c("Acre","Amazonas")
  br_r <- br[br$NAME_1 %in% br_r,]
  br_r <- aggregate(br_r)
  ch <- getData('GADM', country='CHL', level=1)
  ch_r <- c("Arica y Parinacota","Antofagasta","Tarapacá")
  ch_r <- ch[ch$NAME_1 %in% ch_r,]
  ch_r <- aggregate(ch_r)
  ec <- getData('GADM', country='ECU', level=0)
  co <- getData('GADM', country='COL', level=0)
  sa <- bind(pe, br_r, ch_r,ec,co)
  sa <- aggregate(sa)
# #Plot Polygons
  pl.grid <- SpatialPolygons(list(Polygons(list(Polygon(
    cbind(c(-84.47361,-84.47361,-70.37,-70.37), 
          c(-3.3921,-18.41958,-18.41958,-3.3921)),
    FALSE)), '0')), proj4string = CRS(proj4string(sa)))
  poly.only<- gDifference(pl.grid, sa)
  # # #Plot data
plot(poly.only);points(data06t$Lon, data06t$Lat,cex=0.1,col=2)
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
# Project data
poly.onlyutm = spTransform(poly.only, kmproj)
pl.sel = spTransform(pl.grid, kmproj)
## ------------------------------------------------------------------------
mesh <- inla.mesh.2d(boundary = poly.onlyutm,
                     max.edge = c(3,9) * 10,
                         cutoff = 2,# Filter away adjacent points
                         offset = c(30, 150))#
plot(mesh);mesh$n#number of vertices
dim(mesh$graph$vv)#neighborhood
mesh$idx$loc#no one because the mesh is not build according to the points
# water.tri = inla.over_sp_mesh(poly.onlyutm, y = mesh, 
#                               type = "centroid", ignore.CRS = TRUE)#water triangle
# num.tri = length(mesh$graph$tv[, 1])#number of triangles in general
# barrier.tri = setdiff(1:num.tri, water.tri)
# poly.barrier = inla.barrier.polygon(mesh,barrier.triangles = barrier.tri)
x_loc_utm <- spTransform(spdat, kmproj)
#loc.corr <- c(-301.9383,-603.8826 )
#2.Build the spde
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(4, 0.01), # is (range0,Prange) P(range < range0) = Prange
                            prior.sigma = c(1, 0.1)) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
#Prior r0 better not to favour ranges that are smaller than the resolution of the mesh
#Prior sigma how variable is the field from a point to the next
#3. Projector matrix
A= inla.spde.make.A(mesh=mesh, loc=x_loc_utm)
#4.Stack data
stk.dat <- inla.stack(data = list(y = data06t$Anch), 
  A = list(A,1),
  effects = list(list(s = 1:spde$n.spde), 
                 data.frame(Intercept = 1,oceanDist = data06t$DC)),tag = 'dat') 
#5.Formula
f.DC <- y ~ 0 + Intercept + # f is short for formula
  oceanDist +# use the PC prior
  f(s, model = spde)
  
#6.Fitting the model
r.DC <- inla(f.DC, family = 'binomial', # r is short for result
               control.family=list(link='logit'),
               control.compute = list(waic=T,dic=T,cpo=T),
               data = inla.stack.data(stk.dat), 
             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
               control.predictor = list(A = inla.stack.A(stk.dat), link = 1))
summary(r.DC)
#Function scpo to compare models
slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}
#simple  model
onlyDC<-inla(y ~ 0 + Intercept + # f is short for formula
                     oceanDist,family = "binomial",
                     data = inla.stack.data(stk.dat), 
                     control.compute = list(waic=T,dic=T,cpo=T),
                     control.predictor = list(A = inla.stack.A(stk.dat), link = 1))
#compare with spatial
c(only.DC= slcpo(onlyDC), spatialDC = slcpo(r.DC))
#Plot the posterior
par(mfrow = c(2, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.DC$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept', ylab = 'Density') 
plot(r.DC$marginals.fix[[2]], type = 'l', 
     xlab = 'Ocean distance coefficient', ylab = 'Density') 
for (j in 1:2) 
  plot(r.DC$marginals.hy[[j]], type = 'l', 
       ylab = 'Density', xlab = names(r.DC$marginals.hyperpar)[j])
#7.Prediction on a grid
## ----label = "stepsize"--------------------------------------------------
stepsize <- 10 #10km
poly.grid<-spatialEco::remove.holes(poly.only)
poly.gridutm = spTransform(poly.grid, kmproj)# Project data
## ----label = "ncoords"---------------------------------------------------
x.range <- diff(poly.gridutm@bbox[1,])
y.range <- diff(poly.gridutm@bbox[1,])
nxy <- round(c(x.range, y.range) / stepsize)
## ----label = "projgrid"--------------------------------------------------
projgrid <- inla.mesh.projector(mesh, xlim = poly.gridutm@bbox[1,], 
                                ylim = poly.gridutm@bbox[2,], dims = nxy)
## ----label = "projpred"--------------------------------------------------
xmean <- inla.mesh.project(projgrid,r.DC$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, r.DC$summary.random$s$sd)
library(splancs)
xy.in <- inout(projgrid$lattice$loc, as.matrix(poly.gridutm@polygons[[1]]@Polygons[[1]]@coords))#Test points for inclusion in a polygon
table(xy.in)
## ------------------------------------------------------------------------
xmean[!xy.in] <- NA
xsd[!xy.in] <- NA
## ----label = "xrain1", echo = FALSE, fig.width=12, fig.cap = "Posterior mean and standard deviation of the random field (top left and top right, respectively). Posterior mean and standard deviation for the response (bottom left and bottom right, respectively)."----
source('./Scripts/spde-book-functions.R')
#plot of the random field
par(mfrow = c(2, 2), mar = c(0, 0, 0, 4))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xmean))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xsd),
                col = book.color.c2())
#Prediction of the response by computation of posterior distribution
## ------------------------------------------------------------------------
Aprd <- projgrid$proj$A[which(xy.in), ]
##Information to interpolate DC in the grid####
 prdcoo <- projgrid$lattice$loc[which(xy.in), ]#coordinates from the grid onlz in water
# shoreline<-read.csv("F:/Postdoc Hamburg/INLA/shoreline1.csv",sep=",",header = TRUE)
# coordinates(shoreline) <- ~ X + Y
# proj4string(shoreline) <- "+init=epsg:4326"
# shoreline_utm = spTransform(shoreline, kmproj)
# plot(shoreline_utm)
# prdcoodf <-data.frame(prdcoo)
# coordinates(prdcoodf) <- ~ X1 + X2
# proj4string(prdcoodf) <- "+proj=utm +zone=18S ellps=WGS84 +units=km +ellps=WGS84"
# plot(prdcoodf)
#DCgrid <- apply(spDists(shoreline_utm,prdcoodf,longlat = FALSE), 2, min)#only run the first time
#saveRDS(DCgrid, file = "DCgrid1.rds")
DCgrid1 <-readRDS(file = "./Data/DCgrid1.rds")#values for DC for the grid 155*155 but only in water
stk.prd <- inla.stack(
  data = list(y = NA),
  A = list(Aprd, 1), 
  effects = list(s = 1:spde$n.spde, 
                 data.frame(Intercept = 1, oceanDist = DCgrid1)),tag = 'prd') 
stk.all <- inla.stack(stk.dat, stk.prd)#join all in stack
onlyDC2 <- inla(f.DC, family = 'binomial', 
                control.family=list(link='logit'),
                data = inla.stack.data(stk.all), 
                control.predictor = list(A = inla.stack.A(stk.all),
                compute = TRUE, link = 1),
                quantiles = NULL, 
                control.results = list(return.marginals.random = FALSE,
                                       return.marginals.predictor = FALSE), 
                control.mode = list(theta = r.DC$mode$theta,restart = FALSE))
id.prd <- inla.stack.index(stk.all, 'prd')$data#indices of the rows corresponding to the predictions
sd.prd <- m.prd <- matrix(NA, nxy[1], nxy[2])
m.prd[xy.in] <- onlyDC2$summary.fitted.values$mean[id.prd]#mean
sd.prd[xy.in] <- onlyDC2$summary.fitted.values$sd[id.prd]#sd
source('./Scripts/spde-book-functions.R')
par(mfrow = c(2, 2), mar = c(0, 0, 0, 4))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xmean))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xsd),
                col = book.color.c2())
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = m.prd))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = sd.prd),
                col = book.color.c2())
