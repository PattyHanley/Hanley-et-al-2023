require(ncdf4)
require(chron)
require(sp)
require(mgcv)
require(dplyr)
require(cmocean)


BayOfFundy <- nc_open('BayOfFundyBottomTemperature.nc')   # leave this file (pointer) open
latitude <- ncvar_get(BayOfFundy,'lat') 
longitude <- ncvar_get(BayOfFundy,'lon') 
date.str <- as.character(ncvar_get(BayOfFundy,'date_str'))
date <- chron(dates=date.str,format=c(dates="d-mon-y"))
date.num <- as.double(date)
nDays <- length(date.num) - 1  # trimming off empty daily slice -- specific to current 'BayOfFundyBottomTemperature.nc'

basinX <- -66.55
basinY <- 44.40
basin <- c(basinX,basinY)

angle <- 180/180*pi

alngA <- c(1.2*cos(angle),
           1.2*sin(angle))
alngB <- c(-0.1*cos(angle),
           -0.1*sin(angle))
crsA <- c(0.1*cos(angle+pi/2),
          0.1*sin(angle+pi/2))
crsB <- c(-0.1*cos(angle+pi/2),
          -0.1*sin(angle+pi/2))

selected.area.polygon.idx <-
  selected.area.polygon <- 
  rbind( alngA + crsB + basin,
         alngB + crsB + basin,
         alngB + crsA + basin,
         alngA + crsA + basin) %>% 
  round(2)





for (i in 1:dim(selected.area.polygon)[1]) {
  print(i)
  print(c(which(longitude%in%selected.area.polygon[i,1]),
          which(latitude%in%selected.area.polygon[i,2])))
  
selected.area.polygon.idx[i,] <- 
  c(which(longitude%in%selected.area.polygon[i,1]),
    which(latitude%in%selected.area.polygon[i,2]))
}


bathymetry <- ncvar_get(BayOfFundy,'depth')

par(mfrow=c(1,2))

image(longitude,latitude, t(bathymetry),
      xlim=c(-67.2,-66.2),ylim=c(44.25,45.25),asp=c(1,0.8),
      col=topo.colors(100)[seq(100,1,-1)])
contour(longitude,latitude,t(bathymetry),add=TRUE)
polygon(selected.area.polygon,lwd=2)

BayOfFundy.grid.points <- cbind(rep(1:length(longitude),each=length(latitude)),
                                rep(1:length(latitude),length(longitude)))


BayOfFundy.selected.area.idx <- which(in.out(selected.area.polygon.idx, BayOfFundy.grid.points))

# identify depth bin indices to support rapid mapping of temperature to each depth bin
depth.bin.idx <- round(c(bathymetry)[BayOfFundy.selected.area.idx])
depth.bins <- sort(unique(depth.bin.idx))
# print(depth.bins)
depth.bin.indices <- list()
for (b in 1:length(depth.bins)){  # b is for bathymetry, and that's good enough for me
  tmp <-  which(depth.bin.idx == depth.bins[b])
  if(length(tmp) > 0){
    depth.bin.indices[[b]] <- tmp
  }
}  

# set up data arrays
mean.seasonal.temperature.by.depth <- 
  median.seasonal.temperature.by.depth <- 
  sd.seasonal.temperature.by.depth <-
  # min.seasonal.temperature.by.depth <- 
  # max.seasonal.temperature.by.depth <- 
  # Q10.seasonal.temperature.by.depth <- 
  # Q90.seasonal.temperature.by.depth <- 
  array(NA,dim=c(length(depth.bins),nDays))

if(TRUE) {
  for (d in 1:nDays){
    
    print(d)
    temperature <- ncvar_get(BayOfFundy,'temperature',start=c(1,1,d),count=c(-1,-1,1))
    #   print(dim(temperature))
    temperature <- c(temperature)[BayOfFundy.selected.area.idx]
    # contour(temperature)
    for (b in 1:length(depth.bins)){
      # might also want to extract quantile ranges
      temperature.at.depth <- temperature[depth.bin.indices[[b]]]
      if(is.infinite(min(temperature.at.depth,na.rm=TRUE))) print(c(d,b))
      
      # if (length(temperature.at.depth) > 0 ) {
      mean.seasonal.temperature.by.depth[b,d] <- mean(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      median.seasonal.temperature.by.depth[b,d] <- median(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      sd.seasonal.temperature.by.depth[b,d] <- sd(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      # Q10.seasonal.temperature.by.depth[b,d] <- quantile(temperature[depth.bin.indices[[b]]],0.10,na.rm=TRUE)
      # Q90.seasonal.temperature.by.depth[b,d] <- quantile(temperature[depth.bin.indices[[b]]],0.90,na.rm=TRUE)
      # min.seasonal.temperature.by.depth[b,d] <- min(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      # max.seasonal.temperature.by.depth[b,d] <- max(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      # } else {
    }
  }
  
  
  nc_close(BayOfFundy)
  
  # save(file='Central.BayOfFundy.FVCOM.SeasonalBottomTemperature.RData', 
  #      'mean.seasonal.temperature.by.depth','median.seasonal.temperature.by.depth',
  #      # 'min.seasonal.temperature.by.depth','max.seasonal.temperature.by.depth',
  #      # 'Q10.seasonal.temperature.by.depth','Q90.seasonal.temperature.by.depth',
  #      'sd.seasonal.temperature.by.depth',
  #      'depth.bins','date','date.num')
  
  image(1:nDays,depth.bins,t(mean.seasonal.temperature.by.depth),
        col=cmocean('thermal')(100),
        ylim=c(270,0),zlim=c(-3,23),main='mean')
  contour(1:nDays,depth.bins,t(mean.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)

  
}