## Code for model validation - INLA ## 
## Red squat lobster spatiotemporal patterns ## 
## Authors: Noelia Valderrama, Giannina Passuni and Daniel Grados
#################################################################

library(INLA)
library(dplyr)
library(raster)
inla.setOption(pardiso.license = "./Data/pardisovub.lic")
library(maps) # additional helpful mapping packages
source('./Data/spde-book-functions.R')

# Data from sampling 
data = read.csv(file = './Data/datasort_summer_1.csv')
#data<-data%>%filter(Year>2006)
data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)#scaling function
# Standardise covariates
data_1 = data_1 %>% mutate_at(c('DC','SST','SALI','GSST','BAT','CHL'),scale2)
data_1 = na.omit(data_1) #eliminate NAs
data_1$MUN<-ifelse(data_1$MUN>0,1,0)#convert presence absence

# Input data 70 % 30 % 
RNGkind(sample.kind = 'Rounding')
set.seed(1234)
sample_size = floor(0.7*nrow(data_1))
samplet<- sample(seq_len(nrow(data_1)),size = sample_size)
data_1t<-data_1[samplet,]
data_1p<-data_1[-samplet,]
data_1t<-data_1t[order(data_1t$Year),]

#Convert to UTMs
x_loc = cbind(data_1t$Lon_M, data_1t$Lat_M)
spdat <- data_1t
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
x_loc_utm <- spTransform(spdat, kmproj)
x_loc_utm1<-cbind(x_loc_utm@coords[,1],x_loc_utm@coords[,2])

##INLA functions 
#1. General functions for training prediction and validation
#Spatial
mesh = readRDS('./Data/mesh_summer_1.rds')#Mesh from previous steps
mesh$n
spde <- inla.spde2.pcmatern(mesh = mesh,#SPDE set hyperpriors
                            prior.range = c(50, 0.01),
                            prior.sigma = c(5, 0.3))
#Temporal
t<-18
mesh.t1 = inla.mesh.1d(1:18)
h.spec = list(theta=list(prior='pc.cor1', param=c(0, 0.9)))#hyperprior for AR1

#Index space-time
index.st <- inla.spde.make.index('i', n.spde = mesh$n, n.group = 18)

#Formula AR1
f.ar <- z ~ -1 + Intercept + DC + SST + SALI + GSST + BAT + CHL +
  f(i, model = spde,   group = i.group, control.group = list(model = 'ar1', hyper=h.spec))

#Parameters from previous model
#Activate only after training models and used for prediction
#res_ar = readRDS('./Data/r.DCst_spring_ar_0.rds')
#res_ar$mode$theta
#saveRDS(r.ar,'./Data/r.ar.RDS')
r.ar= readRDS('./Data/r.DCst_summer_ar_0.rds')
par_theta<-r.ar$mode$theta
rm(r.ar)
#2. Training#####
#2.1 Projector matrix####
data_1t$Year<-as.numeric(factor(data_1t$Year))
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm1, group.mesh=mesh.t1,
                         group = data_1t$Year)
dim(Aspte)
#2.2 Stack #####
Xm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_1t)
N = nrow(data_1t)
X = data.frame(Intercept = rep(1,N), DC=Xm[,1], SST=Xm[,2], SALI=Xm[,3], 
               GSST=Xm[,4], BAT=Xm[,5],CHL=Xm[,6]) 
#2.2.1 Stack for AR1####
stk.dat.ar <- inla.stack(tag='dat_1', data = list(z = data_1t$MUN),
                         A = list(Aspte,1), effects = list(index.st,X))
#3. Prediction#####
#Data input####
print('Prediction started')
data_pred_1 = read.csv('./Data/var_pred_summer.csv',header = T)
#data_pred_1<-data_pred_1[9166:29328,]
data_pred_1<-data_pred_1%>%dplyr::select(DC,SST,SALI,GSST, BAT,CHL)
# data_pred_1<-data_pred_1%>%mutate(BAT=replace(BAT,BAT>0,NA),CHL=replace(CHL,CHL<0,0))%>%
#   mutate_at(vars(DC:CHL), ~na_if(BAT,NA))#temporal
mean_dat= c(mean(data$DC),mean(data$SST),mean(data$SALI),mean(data$GSST),
            mean(data$BAT),mean(data$CHL))
sd_dat = c(sd(data$DC),sd(data$SST),sd(data$SALI),sd(data$GSST),
           sd(data$BAT),sd(data$CHL))
data_pred = matrix(NaN, ncol = 6, nrow = nrow(data_pred_1))
for (i in 1:6) {
  data_pred[,i] = (data_pred_1[,i]-mean_dat[i])/sd_dat[i]
}
colnames(data_pred) = c('DC','SST','SALI','GSST', 'BAT','CHL')
data_pred = as.data.frame(data_pred)#data scaled
#Coordinates of the mesh
mesh.loc =cbind(rep(mesh$loc[,1],18),rep(mesh$loc[,2],18),sort(rep(seq(1:18),mesh$n)))
#Time for the mesh
time_mesh<-sort(rep(seq(1:18),mesh$n))
#3.1 Projector matrix####
Aprd <- inla.spde.make.A(mesh=mesh,loc=mesh.loc[,1:2],group.mesh=mesh.t1, group = time_mesh)#
dim(Aprd)
#3.2 Stack####
N_pred = nrow(data_pred)
Xp = data.frame(Intercept = rep(1,N_pred), DC=data_pred[,1], SST=data_pred[,2],
                SALI=data_pred[,3], GSST=data_pred[,4], BAT=data_pred[,5],
                CHL=data_pred[,6])
#3.2.1 Stack for AR1####
stk.prd.ar <- inla.stack(tag='pred_1', data = list(z = NA),
                         A = list(Aprd,1), effects = list(index.st,Xp))

#4.Stack all##### 
stk.all.ar <- inla.stack(stk.dat.ar, stk.prd.ar)#join all in stack
#Run INLA#####
r.ar1 <- inla(f.ar,family='binomial', data = inla.stack.data(stk.all.ar),#stk.dat.ar
              control.predictor = list(A = inla.stack.A(stk.all.ar),link=1, compute= TRUE), control.family = list(link='logit'),verbose=TRUE, 
             control.mode = list(theta = par_theta, restart = FALSE),#restart TRUE made crashes INLA for me. Activate control.mode only in prediction
             control.inla=list(strategy='adaptive', int.strategy='eb',verbose=TRUE),control.compute = list(openmp.strategy = 'pardiso.parallel'))#number of threads are the core of the machine you can adjust
#if you are using pardiso  control.compute = list(openmp.strategy = ”pardiso.parallel”)

saveRDS(r.ar1,'pred.ar_summer.rds')
