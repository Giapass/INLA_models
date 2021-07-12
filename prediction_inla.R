library(INLA)
library(dplyr)
library(raster)

# Input data 

data = read.csv(file = './Data/datasort_spring_1.csv')
data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# Standardise covariates
data_1 = data_1 %>% mutate_at(c('DC','SST','SALI','GSST','BAT','CHL'),scale2)

data_1 = na.omit(data_1)

data_1$MUN<-ifelse(data_1$MUN>0,1,0)

# Input data 70 % 30 % 

RNGkind(sample.kind = 'Rounding')
set.seed(1234)
sample_size = floor(0.7*nrow(data_1))
samplet<- sample(seq_len(nrow(data_1)),size = sample_size)
data_1t<-data_1[samplet,]
data_1p<-data_1[-samplet,]

data_1t<-data_1t[order(data_1t$Year),]

x_loc = cbind(data_1t$Lon_M, data_1t$Lat_M)
spdat <- data_1t
coordinates(spdat) <- ~ Lon_M + Lat_M
proj4string(spdat) <- "+init=epsg:4326"
kmproj <- CRS("+proj=utm +zone=18S ellps=WGS84 +units=km")
x_loc_utm <- spTransform(spdat, kmproj)

# Prediction 

print('Prediction started')

mesh = readRDS('./Data/mesh_spring_1.rds')

spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(50, 0.01),
                            prior.sigma = c(5, 0.3))

data_pred = read.csv('./Data/var_pred_spring.csv',header = T)

time_mesh<-sort(rep(seq(1:16),mesh$n))

mesh.loc =cbind(rep(mesh$loc[,1],16),rep(mesh$loc[,2],16),sort(rep(seq(1:16),mesh$n)))

Aprd <- inla.spde.make.A(mesh,loc=mesh.loc,group=time_mesh)#

dim(Aprd)

# #stack prediction

Xpm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_pred)

dim(Xpm)

N_pred = nrow(data_pred)

Xp = data.frame(Intercept = rep(1,N_pred), DC=data_pred[,1], SST=data_pred[,2], SALI=data_pred[,3], GSST=data_pred[,4], BAT=data_pred[,5],CHL=data_pred[,6])

index.st <- inla.spde.make.index('i', n.spde = mesh$n, n.group = 16)

#Replicate
Rep = data_1t$Year

Rep_t = as.data.frame(table(Rep))

Rep_s = rep(c(1:16),times = Rep_t$Freq)

NRep = length(unique(Rep_s))

Arep = inla.spde.make.A(mesh, repl = Rep_s, loc = x_loc_utm)
w.index = inla.spde.make.index(name='w', n.spde = mesh$n, n.repl = NRep)

stk.dat.rp = inla.stack(tag='dat', data = list(y = data_1t$MUN), A = list(Arep, 1), effects = list(w.index,X))

# ar

Xm = model.matrix(~ -1 + DC + SST + SALI + GSST + BAT + CHL, data = data_1t)

N = nrow(data_1t)

X = data.frame(Intercept = rep(1,N), DC=Xm[,1], SST=Xm[,2], SALI=Xm[,3], GSST=Xm[,4], BAT=Xm[,5],CHL=Xm[,6]) 

t<-16
mesh.t1 = inla.mesh.1d(1:16)
data_1t$Year<-as.numeric(factor(data_1t$Year))

Aspte<- inla.spde.make.A(mesh=mesh, loc=x_loc_utm, group.mesh =mesh.t1, group = data_1t$Year)

stk.dat.ar <- inla.stack(tag='dat', data = list(z = data_1t$MUN),
                         A = list(Aspte,1), effects = list(index.st,X))

# stacks 

stk.prd.ar <- inla.stack(tag='pred', data = list(z = NA),
                         A = list(Aprd,1), effects = list(index.st,Xp))

stk.prd.rp <- inla.stack(tag='pred', data = list(z = NA),
                         A = list(Aprd,1), effects = list(w.index,Xp))

stk.all.ar <- inla.stack(stk.dat.ar, stk.prd.ar)#join all in stack

stk.all.rp <- inla.stack(stk.dat.rp, stk.prd.rp)#join all in stack

h.spec = list(theta=list(prior='pc.cor1', param=c(0, 0.9)))

f.ar <- z ~ -1 + Intercept + DC + SST + SALI + GSST + BAT + CHL +f(i, model = spde,   group = i.group, control.group = list(model = 'ar1', hyper=h.spec))

f.rp <- y ~ -1 + Intercept + DC + SST + SALI + GSST + BAT + CHL +f(w, model = spde, replicate =w.repl)

res_ar = readRDS('./Data/r.DCst_spring_ar_0.rds')
res_rp = readRDS('./Data/r.DCst_spring_rp_0.rds')

r.ar <- inla(f.ar,family='binomial', data = inla.stack.data(stk.all.ar), control.predictor = list(A = inla.stack.A(stk.all.ar),link=1, compute= TRUE), control.family = list(link='logit'),verbose=TRUE, control.mode = list(theta = res_ar$mode$theta, restart = FALSE),control.inla=list(strategy='adaptive', int.strategy='eb'))

saveRDS(r.ar,'r.ar_spring.rds')

r.rp <- inla(f.ar,family='binomial', data = inla.stack.data(stk.all.rp), control.predictor = list(A = inla.stack.A(stk.all.rp),link=1, compute= TRUE), control.family = list(link='logit'),verbose=TRUE, control.mode = list(theta = res_rp$mode$theta, restart = FALSE),control.inla=list(strategy='adaptive', int.strategy='eb'))

saveRDS(r.rp,'r.rp_spring.rds')
