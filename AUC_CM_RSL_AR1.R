## Code for AUC and Confusion matrix - INLA ## 
## Red squat lobster spatiotemporal patterns ## 
## Authors: Noelia Valderrama, Giannina Passuni and Daniel Grados
#################################################################

library(INLA)
library(dplyr)
library(raster)


# Data from sampling 

setwd("~/Documents/Master_Tesis/INLA_models-master")

data = read.csv(file = 'datasort_summer_1.csv')
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
data_1p<-data_1p[order(data_1p$Year),]

## second time #######################################

#r.ar= readRDS('./Training/r.DCst_summer_rp_effects_sp.rds')
r.ar1= readRDS('./New_validation/val.ar_summer_DCB.rds')

training = 16651

validation = 7137

ival = seq(training+1,training+validation)

#ival <- inla.stack.index(stk.all.ar,"pred_1")$data
#ival <- inla.stack.index(stk.all.ar,"dat_1")$data

stval <- r.ar1$summary.fitted.values$mean[ival]
#stval <- r.ar$summary.fitted.values$mean[ival]

#summary(data_1p)
#summary(data_1t)

hist(stval)


library(pROC) # install with install.packages("pROC")
library(randomForest) # install with 
#install.packages("randomForest")

roc(data_1p$MUN, stval, plot=TRUE)
##
par(pty = "s") ## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region

## We can calculate the area under the curve...
roc(data_1p$MUN, stval, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", lwd=4, print.auc=TRUE, print.auc.y=40)

stval_1 = rep(1,length(stval))
stval_1[stval<0.5] = 0

#install.packages('caret')
library(caret)

hist(stval)

confusionMatrix(as.factor(stval_1), as.factor(data_1p$MUN))

rm(r.ar1)

