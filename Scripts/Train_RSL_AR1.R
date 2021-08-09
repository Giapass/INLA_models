## Code for model training - INLA ## 
## Red squat lobster spatiotemporal patterns ## 
## Authors: Noelia Valderrama, Giannina Passuni and Daniel Grados
##################################################################

# Libraries 
library (INLA)
library (dplyr)
inla.setOption(pardiso.license = "./Data/pardisovub.lic")
library(raster) # for working with rasters
library(maps) # additional helpful mapping packages
source('./Data/spde-book-functions.R')

# 1. Input data #

data = read.csv(file = './Data/datasort_summer_1.csv')

data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# 2. Standardise covariates
data_1 = data_1 %>% mutate_at(c('DC','SST','SALI','GSST','BAT','CHL'),scale2)

data_1 = na.omit(data_1)

data_1$MUN<-ifelse(data_1$MUN>0,1,0)

# 3. Separating training 70% 

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

# 4. Spatial component - MESH 

mesh = readRDS('./Data/mesh_summer_1.rds')

mesh$n


#4.1 Build the spde

spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(50, 0.01),
                            prior.sigma = c(5, 0.3))


# 5. Spatiotemporal model

t<-18
mesh.t1 = inla.mesh.1d(1:18)
data_1t$Year<-as.numeric(factor(data_1t$Year))

# 5.1 Projector matrix
Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm,group.mesh =mesh.t1, group = data_1t$Year)


# 5.2 Stack data
Xm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_1t)

N = nrow(data_1t)

X = data.frame(Intercept = rep(1,N), DC=Xm[,1], SST=Xm[,2], SALI=Xm[,3], GSST=Xm[,4], BAT=Xm[,5],CHL=Xm[,6])

index.st <- inla.spde.make.index('i', n.spde = mesh$n, n.group = 18)

# 5.3 Index spatiotemporal 
stk.dat.st <- inla.stack(tag='dat', data = list(z = data_1t$MUN),
                         A = list(Aspte,1), effects = list(index.st,X))

# 5.4 Formula

f.DCst <- z ~ -1 + Intercept + DC + SST + SALI + GSST + BAT + CHL +f(i, model = spde,   group = i.group, control.group = list(model = 'ar1', hyper=list(theta=list(prior='pc.cor1', param=c(0, 0.9)))))


# 6. Fitting the model

r.DCst <- inla(f.DCst, family = 'binomial',
               control.family=list(link='logit'),
               data = inla.stack.data(stk.dat.st),
               control.compute = list(waic=TRUE,dic=TRUE,cpo=TRUE,openmp.strategy='huge'),control.predictor = list(A = inla.stack.A(stk.dat.st), link = 1, compute=TRUE), control.inla = list(strategy = 'adaptive', int.strategy = 'eb',tolerance=1e-6),verbose=TRUE)

saveRDS(r.DCst,'r.DCst_summer_ar_0.rds')

