## Code for processing the environmental variables ## 
## Red squat lobster spatiotemporal patterns ## 
## Authors: Noelia Valderrama, Daniel Grados and Giannina Passuni
##################################################################

library(fenix)
library(raster)
library(ncdf4)
library(lubridate)
library(hms)
library(datetime)
library(tidyr)
library(dplyr)
library(stars)
library(ggplot2)
library(sf)


# Data base with presence and absence points

data = read.csv('data.csv')

# 1. Distance to the coast

data$DC = fenix::estima_dc(data$Lon_M, data$Lat_M)

# 2. Bathymetry

# Open the data obtained from GEBCO

bathy = nc_open('gebco_2020_n0.0_s-20.0_w-100.0_e-70.0.nc')

lon = ncvar_get(bathy,"lon")
lat = ncvar_get(bathy,"lat")
bat = ncvar_get(bathy,"elevation")

# Regular grid

Lon = seq(min(lon), max(lon), length.out = length(lon))
Lat = seq(min(lat), max(lat), length.out = length(lat))

XY = merge(Lon,Lat)

XY$bat = matrix(bat,ncol=1)

ras_tot = rasterFromXYZ(XY)

data$BAT = raster::extract(ras_tot, cbind(data$Lon_M, data$Lat_M), method = 'bilinear')

# 3. Classify time (day and night)

Time_M = as.time(data$Time_M)

int <- interval(as_hms(ymd_hms("18:00:00")), as_hms(ymd_hms("05:59:59")))
int2 <- interval(as_hms("06:00:00"), as_hms("17:59:59"))

data$TimeN = Time_M

data$day_night = ifelse(data$TimeN <= as_hms("17:59:59") & data$TimeN >= as_hms("06:00:00"),'Day',ifelse(data$TimeN <= as_hms("05:59:59") | data$TimeN >= as_hms("06:00:00"),'Night',0))


# 4. Classify season (Summer, Spring and cold months)

summer = c('2','3','4')
cold   = c('5','6','8')
spring = c("9" ,"10", "11", '12')

data$season = ifelse(data$Month %in% summer,'Summer',ifelse(data$Month %in% cold,'Cold',ifelse(data_2004$Month %in% spring,'Spring',0)))

# 5. Chlorophyll 

# Open data from Copernicus

nc = nc_open('dataset-oc-glo-bio-multi-l4-chl_interpolated_4km_daily-rep_1614357901760.nc')

# Give a better format to the time

time_chl = ncvar_get(nc,'time')
time_chl = as.Date(time_chl, origin = "1900-01-01")
time_chl = as.character(time_chl)
year_t = substr(time_chl,1,4)
month_t = substr(time_chl,6,7)
day_t = substr(time_chl,9,10)

date_t = paste0(year_t,month_t,day_t)

rm(year_t,month_t,day_t,time_chl)

# Coordinates

lon = ncvar_get(nc, 'lon')
lat = ncvar_get(nc, 'lat')

# Cropping the working area

indlon = which(lon > -86)
indlat = which(lat > -20)

lon = lon[indlon]
lat = lat[indlat]

# Regular grid

Lon = seq(from = min(lon), to = max(lon), length.out = length(lon))
Lat = seq(from = min(lat), to = max(lat), length.out = length(lat))

# Loop for extracting the value in each point of presence and absence

anhos = 2001

for (ianho in anhos)
{
  # Selecting the index corresponding to the year
  
  indanho = which(substr(data$Date_M,1,4) == as.character(ianho))
  
  # Unique dates in the desired year 
  
  Ufechas = unique(data$Date_M[indanho])
  
  # Selecting the value in each point 
  
  for (ifecha in Ufechas)
  {
    print(ifecha)
    
    # Comparison dates from the database (presence and absence) and dates from the environmental database 
    
    indUfecha = which(data$Date_M == ifecha)
    xfecha = which(date_t == ifecha)
    
    ## Chlorophyll from Copernicus has the latitude reverted, to fix this use rev()
    
    chl = ncvar_get(nc,'CHL')[indlon,rev(indlat),xfecha]
  
    # Dataframe and raster
    
    chlXY = merge(Lon,Lat)
    chlXY$z = matrix(chl, ncol = 1)
    chlXY   = rasterFromXYZ(chlXY)
    
    ## Extraction 
    
    result = raster::extract(chlXY, cbind(data$Lon_M[indUfecha],data$Lat_M[indUfecha]), method = 'bilinear', na.rm=T)
    
    data[indUfecha,'chl'] = result
    
  }}

# 6. HYCOM variables ('SST','ELEV','SALI','U','V','GSST') 

# Use DescargaDatos.R and OrdenDownload.R --> Data download

# Verify julian calendar or gregorian

# Extracting variables (Similar code to chlorophyll)

direamb    = '/Volumes/AGAPE/Master_thesis/new_env/'

datos = data

anhos = 2001

# Methodology similar to the one for CHL

for (ianho in anhos)
{
  
  indanho = which( substr(data$Date_M,1,4) == as.character(ianho))
  
  Ufechas = unique(data$Date_M[indanho])
  
  for (ifecha in Ufechas)
  {
    print(ifecha)
    
    indUfecha = which(data$Date_M == ifecha)
    
    carpenc = paste0(direamb,ianho)
    
    nc = nc_open(paste0(carpenc,'/',ifecha,'.nc'))

    lon = ncvar_get(nc, 'lon')
    lat = ncvar_get(nc, 'lat')
    
    # Conversion of coordinates in case the format is different
    
    if (lon[1] > 0) { lon = lon - 360}
    
    indlon = which(lon > -86)
    indlat = which(lat > -20)
    
    lon = lon[indlon]
    lat = lat[indlat]
    
    Lon = seq(from = min(lon), to = max(lon), length.out = length(lon))
    Lat = seq(from = min(lat), to = max(lat), length.out = length(lat))
    
    sst = ncvar_get(nc, 'water_temp')[indlon,indlat,1]
    ele = ncvar_get(nc, 'surf_el')[indlon,indlat]
    sal = ncvar_get(nc, 'salinity')[indlon,indlat,1]
    wvU = ncvar_get(nc, 'water_u')[indlon,indlat,1]
    wvV = ncvar_get(nc, 'water_v')[indlon,indlat,1]
    
    ## Fronts detection 
    
    Gsst = grec::detectFronts(sst)
    
    sstXY = merge(Lon,Lat)
    sstXY$z = matrix(sst, ncol = 1)
    sstXY   = rasterFromXYZ(sstXY)
    
    eleXY = merge(Lon,Lat)
    eleXY$z = matrix(ele, ncol = 1)
    eleXY   = rasterFromXYZ(eleXY)
    
    salXY = merge(Lon,Lat)
    salXY$z = matrix(sal, ncol = 1)
    salXY   = rasterFromXYZ(salXY)
    
    wvUXY = merge(Lon,Lat)
    wvUXY$z = matrix(wvU, ncol = 1)
    wvUXY   = rasterFromXYZ(wvUXY)
    
    wvVXY = merge(Lon,Lat)
    wvVXY$z = matrix(wvV, ncol = 1)
    wvVXY   = rasterFromXYZ(wvVXY)
    
    GsstXY = merge(Lon,Lat)
    GsstXY$z = matrix(Gsst, ncol = 1)
    GsstXY   = rasterFromXYZ(GsstXY)
    
    # Raster Stack
    
    RsStack = stack(sstXY,eleXY,salXY,wvUXY,wvVXY,GsstXY)
    names(RsStack) = c('SST','ELEV','SALI','U','V','GSST')
    
    # Extraction 
  
    result = raster::extract(RsStack, cbind(data$Lon_M[indUfecha],data$Lat_M[indUfecha]), method = 'bilinear')
    
    data[indUfecha,c('SST','ELEV','SALI','U','V','GSST')] = result
    
    rm(RsStack,sstXY,eleXY,salXY,wvUXY,wvVXY,GsstXY,sst,ele,sal,wvU,wvV)
    
  }
  
# Saving per year 
  
    indanhoG = which(substr(data$Date_M,1,4) == ianho)
  
  dataYear = datos[indanhoG,]
  
  nombreguarda = paste0('BaseH_',ianho,'.csv')
  write.csv(dataYear, nombreguarda, row.names = F)
  
  rm(nombreguarda,dataYear)
}

# 7. Direction and speed 

uv2wdws <- function(u,v) {
  
  degrees <- function(radians) 180 * radians / pi
  
  mathdegs <- degrees(atan2(v, u))
  wdcalc <- ifelse(mathdegs > 0, mathdegs, mathdegs + 360)
  wd <- ifelse(wdcalc < 270, 270 - wdcalc, 270 - wdcalc + 360)
  ws <- sqrt(u^2 + v^2)
  
  return(cbind(wd, ws))
  
}

n_wind = uv2wdws(data$U,data$V)
colnames(n_wind) = c('dir','speed')


datos[,c('dir','speed')] = n_wind

