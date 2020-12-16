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
source('./Scripts/spde-book-functions.R')
###Input data#########
an_mu<-read.csv("./Data/data_occ_anc_mun.csv",header=TRUE)
an_mu<-an_mu%>%dplyr::filter(Lon_M>-84.47361)%>%
  mutate(Year=as.factor(format(strptime(Date_M,'%Y%m%d'),'%Y')))%>%drop_na(.)
an_mu<-an_mu%>%dplyr::filter(as.integer(as.character(Year))>=2014)
an_mu$az <- (an_mu$ANC>0) + 0
an_mu$ay<- ifelse(an_mu$ANC>0, an_mu$ANC, NA)
an_mu$mz <- (an_mu$MUN>0) + 0
an_mu$my<- ifelse(an_mu$MUN>0, an_mu$MUN, NA)
#Improve the sampling to have data per year
set.seed(123)
sample_size = floor(0.10*nrow(an_mu))
samplet<- sample(seq_len(nrow(an_mu)),size = sample_size)
an_mut<-an_mu[samplet,]
an_mut<-an_mut[order(an_mut$Year),]
an_mup<-an_mu[-samplet,]
x_loc = cbind(an_mut$Lon_M, an_mut$Lat_M)
spdat <- an_mut
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
##1.Define the boundaries
library(raster)
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
poly.only<-readRDS('F:/Postdoc Hamburg/Master Pleuroncodes/INLA_models/Data/poly.only')
plot(poly.only);points(an_mut$Lon_M, an_mut$Lat_M,cex=0.1,col=2)
# Project data
poly.onlyutm = spTransform(poly.only, kmproj)
## ---------------------Define the mesh---------------------------------------------------
mesh <- inla.mesh.2d(boundary = poly.onlyutm,
                     max.edge = c(4,9) * 10,
                     cutoff = 2,# Filter away adjacent points
                     offset = c(30, 150))#
plot(mesh);mesh$n#number of vertices
dim(mesh$graph$vv)#neighborhood
x_loc_utm <- spTransform(spdat, kmproj)
#2.1 Build the spde
#2.1.1 Spatial priors
prange <- c(20, 0.9)#
psigma<-c(4, 0.5)#
spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = prange, # is (range0,Prange) P(range < range0) = Prange
                            prior.sigma = psigma) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
#2.2 Define the time effects
t<-6
mesh.t1 = inla.mesh.1d(1:t)
time<-as.numeric(factor(an_mut$Year))
#3. Projector matrix#######
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm,
                         group.mesh = mesh.t1, group = time)#A dim[mesh,location*time]
dim(Aspte)==c(nrow(an_mut),(t*spde$n.spde))#to confirm the correct dim of matrix A
#4.Stack data##############
#4.1 Make index NEW
rf1 <- inla.spde.make.index('rf1', n.spde = spde$n.spde,
                                  n.group = t)#INDEX SPATIO TEMPORAL FROM BINOMIAL ANCHOVY
rf2 <- inla.spde.make.index('rf2', n.spde = spde$n.spde,
                            n.group = t)#INDEX SPATIO TEMPORAL FROM BINOMIAL MUNIDA
brf1 <- inla.spde.make.index('brf1', n.spde = spde$n.spde,
                                  n.group = t)#NEW INDEX SPATIO TEMPORAL SHARED INTRASPECIFIC ANCHOVY
brf2 <- inla.spde.make.index('brf2', n.spde = spde$n.spde,
                             n.group = t)#INDEX SPATIO TEMPORAL SHARED INTRASPECIFIC MUNIDA
g1rf1 <- inla.spde.make.index('g1rf1', n.spde = spde$n.spde,
                            n.group = t)#INDEX SPATIO TEMPORAL FROM BINOMIAL INTERSPECIES
g2rf1 <- inla.spde.make.index('g2rf1', n.spde = spde$n.spde,
                             n.group = t)#INDEX SPATIO TEMPORAL FROM GAMMA INTERSPECIES
#4.2 Stack estimation
stk.zan<- inla.stack(data = list(Y = cbind(an_mut$az,NA,NA,NA),link=1),#NEW 
                        A = list(Aspte),
                        effects = list(rf1=rf1),tag = 'an.z') 
stk.yan <- inla.stack(data = list(Y = cbind(NA,an_mut$ay,NA,NA),link=2),#NEW 
                        A = list(Aspte),
                        effects = list(brf1),tag = 'an.y') 
stk.zmu <- inla.stack(data=list(Y = cbind(NA,NA,an_mut$mz,NA),link=3),
                      A=list(Aspte,Aspte,1), tag='mu.z',
                      effects=list(g1rf1 = g1rf1,rf2 = rf2,
                                   alpha1=rep(1,length(an_mut$mz))))
stk.ymu<-inla.stack(data=list(Y = cbind(NA,NA,NA,an_mut$my),link=4),tag="mu.y",
                    A=list(Aspte,Aspte,1),
                    effects=list(g2rf1 = g2rf1,brf2 = brf2,
                                 alpha2=rep(1,length(an_mut$my))))
stack<-inla.stack(stk.zan,stk.yan,stk.zmu,stk.ymu)
#Hyperparameters
#intra specific
beta_1 = list(theta=list(prior='normal', param=c(0,1)))
beta_2 = list(theta=list(prior='normal', param=c(0,1)))
#inter specific
beta_z2 = list(theta=list(prior='normal', param=c(0,1)))
gamma_z2 = list(theta=list(prior='normal', param=c(0,1)))
#####Priors######################
rhoprior <- list(theta = list(prior = 'pccor1',param = c(0.5, 0.9)))#rho
pcgprior <- list(prior = 'pc.gamma', param = 1)#gamma prior
#5.Formula New
f.interaction <- Y ~ -1+ alpha1+alpha2+
  f(rf1,model=spde, group=rf1.group, control.group = list(model = 'ar1', hyper = rhoprior))+#binomial anchovy
  f(brf1, copy = "rf1",fixed=FALSE,group=brf1.group,hyper = list(beta_1))+#gamma anchovy
  f(g1rf1,copy = "rf1",fixed=FALSE,group=g1rf1.group,hyper = list(beta_z2))+#shared rf binomial
  f(g2rf1,copy = "rf1",fixed=FALSE,group=g2rf1.group, hyper = list(gamma_z2))+#shared rf gamma
  f(rf2, model = spde,group = rf2.group,control.group = list(model = 'ar1', hyper = rhoprior)) +#binomial munida
  f(brf2, copy = 'rf2' , fixed = FALSE,group=brf2.group,hyper = list(beta_2))#gamma munida
#6.Fitting the model
r.interaction <- inla(f.interaction, data = inla.stack.data(stack), 
                  family =rep(c("binomial", "gamma"),2), 
                  control.predictor = list(A = inla.stack.A(stack), link = 1,compute=TRUE),
                  control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),num.threads=2,
                  verbose=TRUE,control.compute = list(dic=TRUE))
summary(r.interaction )



#Hurdle model munida
#2.1 Build the spde
#2.1.1 Spatial priors
prange <- c(20, 0.9)#
psigma<-c(4, 0.5)#
rhoprior <- list(theta = list(prior = 'pccor1',param = c(0.5, 0.9)))#rho
pcgprior <- list(prior = 'pc.gamma', param = 1)#gamma prior
spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = prange, # is (range0,Prange) P(range < range0) = Prange
                            prior.sigma = psigma) # is (sigma0,Psigma) P(sigma > sigma_0) = Psigma
#2.2 Define the time effects
t<-6
mesh.t1 = inla.mesh.1d(1:t)
time<-as.numeric(factor(an_mut$Year))
#3. Projector matrix#######
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm,
                         group.mesh = mesh.t1, group = time)#A dim[mesh,location*time]
dim(Aspte)==c(nrow(an_mut),(t*spde$n.spde))#to confirm the correct dim of matrix A
#4.Stack data##############
#4.1 Make index NEW
rf2 <- inla.spde.make.index('rf2', n.spde = spde$n.spde,
                            n.group = t)#INDEX SPATIO TEMPORAL FROM BINOMIAL MUNIDA
brf2 <- inla.spde.make.index('brf2', n.spde = spde$n.spde,
                             n.group = t)#INDEX SPATIO TEMPORAL SHARED INTRASPECIFIC MUNIDA

#Hyperparameters
#intra specific
beta_2 = list(theta=list(prior='normal', param=c(0,1)))

stk.zmu <- inla.stack(data=list(Y = cbind(an_mut$mz,NA),link=1),
                      A=list(Aspte,1), tag='mu.z',
                      effects=list(rf2,
                                   alpha1=rep(1,length(an_mut$mz))))
stk.ymu<-inla.stack(data=list(Y = cbind(NA,an_mut$my),link=2),tag="mu.y",
                    A=list(Aspte,1),
                    effects=list(brf2,
                                 alpha2=rep(1,length(an_mut$my))))
stack<-inla.stack(stk.zmu,stk.ymu)
f.interaction <- Y ~ -1+ alpha1+alpha2+
  f(rf2, model = spde,group = rf2.group,control.group = list(model = 'ar1', hyper = rhoprior)) +#binomial munida
  f(brf2, copy = 'rf2' , fixed = FALSE,group=brf2.group,hyper = list(beta_2))#gamma munida
#6.Fitting the model
r.interaction <- inla(f.interaction, data = inla.stack.data(stack), 
                      family =c("binomial", "gamma"), 
                      control.predictor = list(A = inla.stack.A(stack), link = 1,compute=TRUE),
                      control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),num.threads=2,
                      verbose=TRUE,control.compute = list(dic=TRUE))
#Exlplore the marginals
par(mfrow = c(4, 2), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
plot(r.interaction$marginals.fix[[1]], type = 'l', 
     xlab = 'Intercept.Binomial', ylab = 'Density')
plot(r.interaction$marginals.fix[[2]], type = 'l', 
     xlab = 'Intercept.Gamma', ylab = 'Density')
par(mfrow = c(3, 3), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 
for (j in 1:5) 
  plot(r.interaction$marginals.hy[[j]], type = 'l', xlab = names(r.interaction$marginals.hyperpar)[j])
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
  r <- inla.mesh.project(projgrid,field = r.interaction$summary.ran$rf2$mean[idx])
  r[!xy.in] <- NA
  return(r)
})

mu.st2 <- lapply(1:mesh$n, function(j) {
  idx <- 1:spde$n.spde + (j - 1) * spde$n.spde
  r <- inla.mesh.project(projgrid,field = r.interaction$summary.ran$brf2$mean[idx])
  r[!xy.in] <- NA
  return(r)
})
par(mfrow = c(3, 2), mar = c(0, 0, 1, 0))
zlm1 <- range(unlist(mu.st1), na.rm = TRUE)
# identify wich time location is near each knot
for (j in 1:6) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st1[[j]]), 
                  zlim = zlm1, main = paste0("Mesh timepoint: ", j))
}
for (j in 1:6) {
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = mu.st2[[j]]), 
                  main = paste0("Mesh timepoint: ", j))
}

