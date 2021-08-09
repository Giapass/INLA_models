## Code for VIF - INLA ## 
## Red squat lobster spatiotemporal patterns ## 
## Authors: Noelia Valderrama, Giannina Passuni and Daniel Grados
#################################################################

setwd("~/Documents/Master_Tesis/INLA_models-master")


library(tidyverse)
library(caret)


# Load the data
data = read.csv(file = 'datasort_summer_1.csv')

data_1<-data%>%dplyr::select(Year,Lon_M,Lat_M,MUN,ANC,DC,SST,SALI,GSST, BAT,CHL)
# Standardise covariates
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

# Build the model
model1 <- lm(MUN ~DC+SST+SALI+GSST+BAT+CHL, data = data_1t)
# Make predictions
predictions <- model1 %>% predict(data_1p)
# Model performance
data.frame(
  RMSE = RMSE(predictions, data_1p$MUN),
  R2 = R2(predictions, data_1p$MUN))

car::vif(model1)

