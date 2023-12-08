require(ncdf4)
require(chron)
require(sp)
require(mgcv)


BayOfFundy <- nc_open('BayOfFundyBottomTemperature.nc')   # leave this file (pointer) open
latitude <- ncvar_get(BayOfFundy,'lat') 
longitude <- ncvar_get(BayOfFundy,'lon') 
date.str <- as.character(ncvar_get(BayOfFundy,'date_str'))
date <- chron(dates=date.str,format=c(dates="d-mon-y"))
date.num <- as.double(date)
nDays <- length(date.num) - 1  # trimming off empty daily slice -- specific to current 'BayOfFundyBottomTemperature.nc'



# read in bottom temperature field for each day
# calculate statistics for temperatures observed in each depth bin


# define lon-lat points of an arbitrary polygon, ensuring progress "around a circle" so there are no "twists" in the resulting shape
# i.e., crossing edges;  values here are for the polygon we originally used, shifted a bit to the east

# selected.area.polygon.idx <- 
#   selected.area.polygon <- rbind( c(-65.0,44.0),  # Original square shifted east to better match inferred movement ranges of tagged lobsters
#                                 c(-65.0,45.5),    # Note that you could do this with different areas to see whether there are possible
#                                 c(-67.5,45.5),    # differences in lobster behavior (movement tracks) conditional on location
#                                 c(-67.5,44.0))    # for the same movement rules


selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-64.5,45.0),  # Possible quadrilateral to section bay across its width
                                  c(-65.5,43.75),    
                                  c(-68.0,44.75),    
                                  c(-66.0,46.0))    

selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-65.85,45.0),  # Possible quadrilateral to section bay across its width
                                  c(-66.65,44.35),    
                                  c(-68.0,44.75),    
                                  c(-66.5,45.5)) 

selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-65.85,45.05),  # Possible quadrilateral to section bay across its width
                                  c(-66.65,44.35),    
                                  c(-67.25,44.75),    
                                  c(-66.0,45.35)) 

selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-65.85,45.05),  # Possible quadrilateral to section bay across its width
                                  c(-66.65,44.35),    
                                  c(-66.00,44.00),    
                                  c(-65.55,44.75)) 
selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-64.5,45.0),  
                                  c(-65.5,43.5),    
                                  c(-68.4,44.5),    
                                  c(-66.0,45.75) )  

selected.area.polygon.idx <-
  selected.area.polygon <- rbind(c(-66.45, 44.3),
                                 c(-66.7, 44.35),
                                 c(-66.8, 44.5),
                                 c(-66.79,44.82),
                                 c(-66.37, 44.7))

selected.area.polygon.idx <- 
  selected.area.polygon <- rbind( c(-66.8,44.8),  # Possible quadrilateral to section bay across its width
                                  c(-66.8, 44.58),    
                                  c(-66.3, 44.58),
                                  c(-66.3, 44.8)
  )  


for (i in 1:dim(selected.area.polygon)[1]) {
  print(i)
  print(c(which(longitude%in%selected.area.polygon[i,1]),which(latitude%in%selected.area.polygon[i,2])))
  selected.area.polygon.idx[i,] <- c(which(longitude%in%selected.area.polygon[i,1]),which(latitude%in%selected.area.polygon[i,2]))
}


bathymetry <- ncvar_get(BayOfFundy,'depth')

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
  min.seasonal.temperature.by.depth <- 
  max.seasonal.temperature.by.depth <- 
  Q10.seasonal.temperature.by.depth <- 
  Q90.seasonal.temperature.by.depth <- 
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
      Q10.seasonal.temperature.by.depth[b,d] <- quantile(temperature[depth.bin.indices[[b]]],0.10,na.rm=TRUE)
      Q90.seasonal.temperature.by.depth[b,d] <- quantile(temperature[depth.bin.indices[[b]]],0.90,na.rm=TRUE)
      min.seasonal.temperature.by.depth[b,d] <- min(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      max.seasonal.temperature.by.depth[b,d] <- max(temperature[depth.bin.indices[[b]]],na.rm=TRUE)
      # } else {
    }
  }
  
  
  nc_close(BayOfFundy)
  
  save(file='Central.BayOfFundy.FVCOM.SeasonalBottomTemperature.RData', 
       'mean.seasonal.temperature.by.depth','median.seasonal.temperature.by.depth',
       'min.seasonal.temperature.by.depth','max.seasonal.temperature.by.depth',
       'Q10.seasonal.temperature.by.depth','Q90.seasonal.temperature.by.depth',
       'sd.seasonal.temperature.by.depth',
       'depth.bins','date','date.num')
  
  image(1:nDays,depth.bins,t(mean.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='mean')
  contour(1:nDays,depth.bins,t(mean.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  # image(1:nDays,depth.bins,t(median.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='median')
  # contour(1:nDays,depth.bins,t(median.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  image(1:nDays,depth.bins,t(min.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='minimum')
  contour(1:nDays,depth.bins,t(min.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  # image(1:nDays,depth.bins,t(max.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='maximum')
  # contour(1:nDays,depth.bins,t(max.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  # image(1:nDays,depth.bins,t(sd.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(0,4),main='standard deviation')
  # contour(1:nDays,depth.bins,t(sd.seasonal.temperature.by.depth),levels=seq(0.5,4,0.5),add=TRUE)
  # 
  # image(1:nDays,depth.bins,t(Q10.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='quantile 10')
  # contour(1:nDays,depth.bins,t(Q10.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  # image(1:nDays,depth.bins,t(Q90.seasonal.temperature.by.depth),col=rainbow(100),ylim=c(270,0),zlim=c(-3,23),main='quantile 90')
  # contour(1:nDays,depth.bins,t(Q90.seasonal.temperature.by.depth),levels=seq(0,20,2),add=TRUE)
  
}