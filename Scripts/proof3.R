
library (INLA)
INLA:::inla.dynload.workaround()
library (dplyr)
library(rgeos)
library(raster) # for working with rasters
library(maps) # additional helpful mapping packages
library(maptools)
source('./Scripts/spde-book-functions.R')

# Input data #

# Logaritmizar CHL 

data = read.csv(file = './Data/data_filtered_complete.csv')

data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)
data_1 = na.omit(data_1)
data_1$MUN<-ifelse(data_1$MUN>0,1,0)
data_1$BAT<-(data_1$BAT*-1)/1000
# Standardise covariates, NOT NEEDED, doesnt help to optmise computing time
#data_1<-data_1%>%mutate_at(vars(DC,SST,SALI,GSST, BAT,CHL),scale)


RNGkind(sample.kind = 'Rounding')
set.seed(1234)
sample_size = floor(0.7*nrow(data_1))
samplet<- sample(seq_len(nrow(data_1)),size = sample_size)
data_1t<-data_1[samplet,]
data_1p<-data_1[-samplet,]

data_1t<-data_1t[order(data_1t$Year),]#important to reorder
table(data_1t$Year)
x_loc = cbind(data_1t$Lon_M, data_1t$Lat_M)
spdat <- data_1t
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")

# Project data
x_loc_utm <- spTransform(spdat, kmproj)

# mesh

mesh = readRDS('./Data/mesh.rds') 

mesh$n

#2.1 Build the spde
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(60, 0.1), 
                            prior.sigma = c(1, 0.01)) 
 

t<-18
mesh.t1 = inla.mesh.1d(1:18)
data_1t$Year<-as.numeric(factor(data_1t$Year))

#3. Projector matrix
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm, group.mesh=mesh.t1,group = data_1t$Year)

dim(Aspte)

#4.Stack data
Xm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_1t)

N = nrow(data_1t)

X = data.frame(Intercept = rep(1,N), DC=Xm[,1], SST=Xm[,2], SALI=Xm[,3], GSST=Xm[,4], BAT=Xm[,5],CHL=Xm[,6]) 

index.st <- inla.spde.make.index('i', n.spde = mesh$n, n.group = 18)

#NEW INDEX SPATIO TEMPORAL
stk.dat.st <- inla.stack(tag='dat', data = list(z = data_1t$MUN), 
                         A = list(Aspte,1), effects = list(index.st,X)) 
#Another way to noted the fixed effects
#stk.dat.st <- inla.stack(tag='dat', data = list(z = data_1t$MUN), 
#                         A = list(Aspte,1), effects = list(index.st,list(Intercept=rep(1,N),
#                         DC=data_1t$DC,SST=data_1t$SST,SALI=data_1t$SALI,GSST=data_1t$GSST,CHL=data_1t$CHL))) 

# #4.1 stack for prediction
# #A projector for prediction
time_mesh<-sort(rep(seq(1:18),mesh$n))

mesh.loc =cbind(rep(mesh$loc[,1],18),rep(mesh$loc[,2],18),sort(rep(seq(1:18),mesh$n)))

Aprd <- inla.spde.make.A(mesh,loc=mesh.loc,n.group=mesh.t1$n,group=time_mesh)#

dim(Aprd)

data_pred = read.csv('./Data/var_pred_c.csv',header = T)


#Stack data
Xpm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_pred)

dim(Xpm)

N_pred = nrow(data_pred)

Xp = data.frame(Intercept = rep(1,N_pred), DC=data_pred[,1], SST=data_pred[,2], SALI=data_pred[,3], GSST=data_pred[,4], BAT=data_pred[,5],CHL=data_pred[,6]) 


# # #Mesh location
# #stack prediction
stk.prd.st <- inla.stack(tag='pred', data = list(z = NA),
                         A = list(Aprd,1), effects = list(index.st,Xp))
#Another way to noted the fixed effects
#stk.prd.st <- inla.stack(tag='pred', data = list(z = NA),
#                         A = list(Aprd,1), effects = list(index.st,list(Intercept=rep(1,N),
#                         DC=data_pred$DC,SST=data_pred$SST,SALI=data_pred$SALI,GSST=data_pred$GSST,CHL=data_pred$CHL)))
#4.2 stack all
stk.all.st <- inla.stack(stk.dat.st,stk.prd.st)

#5.Formula
f.DCst <- z ~ -1 + Intercept + DC + SST + SALI + GSST + BAT + CHL +f(i, model = spde,   group = i.group, control.group = list(model = 'ar1',
hyper=list(theta=list(prior='pc.cor1', param=c(0.3, 0.7)))))

#6.Fitting the model
r.DCst <- inla(f.DCst, family = 'binomial', 
               control.family=list(link='logit'),
               data = inla.stack.data(stk.all.st), 
               control.predictor = list(A = inla.stack.A(stk.all.st), link = 1, compute=TRUE), control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),verbose=TRUE)

#summary(r.DCst)

#saveRDS(r.DCst,'proof3.rds')

r.DCst = readRDS('proof3.rds')

summary(r.DCst)


#Plot the posterior

par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.DCst$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept', ylab = 'Density') 
plot(r.DCst$marginals.fix[[2]], type = 'l', 
     xlab = 'Ocean distance coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[3]], type = 'l', 
     xlab = 'SST coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[4]], type = 'l', 
     xlab = 'Salinity coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[5]], type = 'l', 
     xlab = 'GSST coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[6]], type = 'l', 
     xlab = 'Bathymetry coefficient', ylab = 'Density') 
plot(r.DCst$marginals.fix[[7]], type = 'l', 
     xlab = 'CHL coefficient', ylab = 'Density') 


par(mfrow = c(3, 1), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 

for (j in 1:3) 
  plot(r.DCst$marginals.hy[[j]], type = 'l', 
       ylab = 'Density', xlab = names(r.DCst$marginals.hyperpar)[j])
#7.Prediction on a grid
## ----GRID for projection of field --------------------------------------------------
#GRID 10km###
stepsize <- 10 #10km
poly.only = readRDS('Data/poly.only')
poly.grid<-spatialEco::remove.holes(poly.only)
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
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
  idx <- 1:spde$n.spde + (j - 1) * spde$n.spde #index to find the spatial field per year
  r <- inla.mesh.project(projgrid,field = r.DCst$summary.random$i$mean[idx])
  r[!xy.in] <- NA
  return(r)
})
zlm1 <- range(unlist(mu.st1), na.rm = TRUE)
# plot of the random field
x11()
par(mfrow = c(6, 3), mar = c(0, 0, 1, 0))

for (j in 1:18) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st1[[j]]), 
                  zlim = zlm1, main = paste0("Mesh timepoint: ", j))
}

#Prediction of the response by computation of posterior distribution
## ------------------------------------------------------------------------
source('./Scripts/spde-book-functions.R')
# #Matrix for prediction
index_predic <- inla.stack.index(stk.all.st,"prd")$data
stpred <- matrix(r.DCst$summary.fitted.values$mean[index_predic], 
                 spde$n.spde)
x11()
par(mfrow = c(6, 3), mar = c(1, 1, 1, 1))
for(i in 1:18){
  meanproj <- inla.mesh.project(projgrid,field=stpred[,i])
  meanproj[!xy.in] <- NA
  print(mean(meanproj,na.rm=TRUE))
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = meanproj))
}

#Plot per year spatial field

