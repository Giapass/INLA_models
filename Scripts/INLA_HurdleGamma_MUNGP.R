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
data_1$z <- (data_1$MUN>0) + 0
data_1$y<- ifelse(data_1$MUN>0,data_1$MUN, NA)
data_1$BAT<-(data_1$BAT*-1)/1000
RNGkind(sample.kind = 'Rounding')
set.seed(1234)
sample_size = floor(0.7*nrow(data_1))
samplet<- sample(seq_len(nrow(data_1)),size = sample_size)
data_1t<-data_1[samplet,]
data_1p<-data_1[-samplet,]

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

#####Priors######################
prange <- c(60, 0.1)#spde
psigma<-c(1, 0.01)#spde
rhoprior <- list(theta = list(prior = 'pccor1',param = c(0.3, 0.7)))#rho
bprior <- list(prior = 'gaussian', param = c(0,1))#betaprior, shared
pcgprior <- list(prior = 'pc.gamma', param = 1)#gamma prior

#2.1 Build the spde
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = prange, # is (range0,Prange) P(range < range0) = Prange ## 25 km instead of 40 
                            prior.sigma = psigma) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
# 
#Prior r0 better not to favour ranges that are smaller than the resolution of the mesh
#Prior sigma how variable is the field from a point to the next

#2.2 Define the time effects,NEW
t<-18#18
mesh.t1 = inla.mesh.1d(1:18)
data_1t$Year<-as.numeric(factor(data_1t$Year))
#3. Projector matrix
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm, group.mesh = mesh.t1,#anadi group.mesh
                         group = data_1t$Year)#NEW
dim(Aspte)
#4.Stack data##############
#4.1 Make index
index.stz <- inla.spde.make.index('i_z', n.spde = spde$n.spde,
                                  n.group = t)#NEW INDEX SPATIO TEMPORAL
index.stc <- inla.spde.make.index('i_c', n.spde = spde$n.spde,
                                  n.group = t)#NEW INDEX SPATIO TEMPORAL
index.sty <- inla.spde.make.index('i_y', n.spde = spde$n.spde,
                                  n.group = t)#NEW INDEX SPATIO TEMPORAL
#INDEX SPATIO TEMPORAL
N = nrow(data_1t)
#4.2 Stack estimation
stk.z.est <- inla.stack(data = list(Y = cbind(data_1t$z,NA),link=1),#NEW 
                        A = list(Aspte,1),
                        effects = list(index.stz, 
                                       data.frame(z.Intercept=rep(1,N),
                                                  DC_z= data_1t$DC,BAT_z=data_1t$BAT,SST_z=data_1t$SST,
                                                  SALI_z=data_1t$SALI,CHL_z=data_1t$CHL,GSST_z=data_1t$GSST)),
                        tag = 'zdat') 
stk.y.est <- inla.stack(data = list(Y = cbind(NA,data_1t$y),link=2),#NEW 
                        A = list(Aspte,1),
                        effects = list(c(index.stc,index.sty), 
                                       data.frame(y.Intercept = rep(1,N),
                                                  DC_y = data_1t$DC,BAT_y=data_1t$BAT,SST_y=data_1t$SST,
                                                  SALI_y=data_1t$SALI,CHL_y=data_1t$CHL,GSST_y=data_1t$GSST)),
                        tag = 'ydat') 
#4.3. Predictions
time_mesh<-sort(rep(seq(1:18),mesh$n))
mesh.loc =cbind(rep(mesh$loc[,1],18),rep(mesh$loc[,2],18),sort(rep(seq(1:18),mesh$n)))
meshcoord <-data.frame(mesh$loc[,1],mesh$loc[,2])
names(meshcoord)<-c('utmx','utmy')
coordinates(meshcoord) <- ~ utmx + utmy
proj4string(meshcoord) <- "+proj=utm +zone=18S ellps=WGS84 +units=km +ellps=WGS84"
Aprd <- inla.spde.make.A(mesh,loc=mesh.loc,n.group=mesh.t1$n,group=time_mesh)#
dim(Aprd)
data_pred = read.csv('var_pred_c.csv')
#Variables
#Xpm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_pred)
N_pred = nrow(data_pred)
#Xp = data.frame(Intercept = rep(1,N_pred), DC=Xpm[,1], SST=Xpm[,2], SALI=Xpm[,3], GSST=Xpm[,4], BAT=Xpm[,5],CHL=Xpm[,6]) 
stk.z.prd <- inla.stack(data = list(Y = matrix(NA, ncol(Aprd), 2), link = 1),
                        A = list(Aprd,1), 
                        effects = list(index.stz,
                                       data.frame(z.Intercept=rep(1,N_pred),
                                                  DC_z= data_pred$DC,BAT_z=data_pred$BAT,SST_z=data_pred$SST,
                                                  SALI_z=data_pred$SALI,CHL_z=data_pred$CHL,GSST_z=data_pred$GSST)),tag = 'zprd') 
stk.y.prd <- inla.stack(data = list(Y = matrix(NA, ncol(Aprd), 2), link = 2),
                        A = list(Aprd,1), 
                        effects = list(c(index.stc,index.sty),data.frame(y.Intercept = rep(1,N_pred),
                                                                         DC_y = data_pred$DC,BAT_y=data_pred$BAT,SST_y=data_pred$SST,
                                                                         SALI_y=data_pred$SALI,CHL_y=data_pred$CHL,GSST_y=data_pred$GSST)),tag = 'yprd') 

#4.4 Stack all
stk.all.hurdle <- inla.stack(stk.z.est,stk.y.est, stk.z.prd, stk.y.prd)#

#5.Formula 
#5.1 Shared space time
f.DCjoint <- Y ~ -1 + z.Intercept +y.Intercept+DC_z+DC_y+BAT_z+BAT_y+SST_z+SST_y+SALI_z+SALI_y+CHL_z+CHL_y+GSST_z+GSST_y+
  f(i_z,model=spde, group=i_z.group, control.group = list(model = 'ar1', hyper = rhoprior))+
  f(i_c, copy = "i_z", fixed=FALSE,group = i_c.group, hyper = list(theta = bprior)) + #fixed false to have a parameter of the shared space-time
  f(i_y, model = spde, group = i_y.group, control.group =list(model = 'ar1', hyper = rhoprior)) # bprior is a term of interaction between presence and density 

#6.Fitting the model
#Define link as a vector with information for every point in estimation and prediction
#linkt<-c(rep(1,nrow(data_1t)),rep(2,nrow(data_1t)))#Noelia, no lo necesitas es solo mi ejercicio de hacer nuevas estructuras
#linkt<-c(rep(1,nrow(data_1t)),rep(2,nrow(data_1t)),rep(1,mesh$n*t)),rep(2,mesh$n*t))#Noelia, no lo necesitas es solo mi ejercicio de hacer nuevas estructuras
r.DCjoint <- inla(f.DCjoint, data = inla.stack.data(stk.all.hurdle), 
                  family =c("binomial", "gamma"), 
                  control.family =list(list(), list(hyper = list(theta = pcgprior))), 
                  control.predictor = list(A = inla.stack.A(stk.all.hurdle), 
                                           link = inla.stack.data(stk.all.hurdle)$link,compute=TRUE),
                  control.inla = list(strategy = 'adaptive', int.strategy = 'eb',cmin=0),verbose=TRUE)
summary(r.DCjoint)
saveRDS(r.DCjoint,'hurdle_shared_1119.rds')
#Estimation
r.DCjoint$summary.fitted.values[inla.stack.index(stk.all.hurdle,'ydat')$data,'mean']#save them before going out of the cluster
r.DCjoint$summary.fitted.values[inla.stack.index(stk.all.hurdle,'zdat')$data,'mean']#save them before going out of the cluster
dat4bp_est<-data.frame(group=sort(rep(seq(1:5),nrow(data_1t))),
                       ppres=r.DCjoint$summary.fitted.values[inla.stack.index(stk.all.hurdle,'zdat')$data,'0.5quant'])
boxplot(ppres~group,data=dat4bp_est)#TERRIBLE, it means we should use prediction to obtain the maps
#Prediction
#r.DCjoint$summary.fitted.values[inla.stack.index(stk.all.hurdle,'ypred')$data,'mean']#save them before going out of the cluster
#r.DCjoint$summary.fitted.values[inla.stack.index(stk.all.hurdle,'zpred')$data,'mean']#save them before going out of the cluster
#Plot the posterior
plot(inv.logit(r.DCjoint$summary.random$i_z$`0.5quant`))
plot(exp(r.DCjoint$summary.random$i_y$`0.5quant`))
plot(exp(r.DCjoint$summary.random$i_c$`0.5quant`))
boxplot(ppres~group,data=dat4bp)
dat4bp<-data.frame(group=sort(rep(seq(1:18),3775)),ppres=inv.logit(r.DCjoint$summary.random$i_z$mean))
boxplot(dat4bp$group)
par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.DCjoint$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept.Z', ylab = 'Density')
plot(r.DCjoint$marginals.fix[[2]], type = 'l', 
     xlab = 'Intercept.y', ylab = 'Density')
plot(r.DCjoint$marginals.fix[[3]], type = 'l', 
     xlab = 'Ocean distance coefficient_z', ylab = 'Density') 
plot(r.DCjoint$marginals.fix[[4]], type = 'l', 
     xlab = 'Ocean distance coefficient_y', ylab = 'Density') 
par(mfrow = c(3, 3), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
for (j in 1:7) 
  plot(r.DCjoint$marginals.hy[[j]], type = 'l', xlab = names(r.DCjoint$marginals.hyperpar)[j])
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
  r <- inla.mesh.project(projgrid,field = r.DCjoint$summary.ran$i_z$mean[idx])
  r[!xy.in] <- NA
  return(r)
})

mu.st2 <- lapply(1:mesh$n, function(j) {
  idx <- 1:spde$n.spde + (j - 1) * spde$n.spde
  r <- inla.mesh.project(projgrid,field = r.DCjoint$summary.ran$i_y$mean[idx])
  r[!xy.in] <- NA
  return(r)
})
par(mfrow = c(3, 2), mar = c(0, 0, 1, 0))
zlm1 <- range(unlist(mu.st1), na.rm = TRUE)
# identify which time location is near each knot
for (j in 1:18) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st1[[j]]), 
                  zlim = zlm1, main = paste0("Binomial: ", years[j]))
}
for (j in 1:18) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st2[[j]]), 
                  main = paste0("Gamma: ",years[j]))
}