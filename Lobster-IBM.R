Lobster.1D.IBM <- function(n.lobsters = 1000,          # number of lobster in simulation
                           useTemperatureRule = TRUE, # TRUE = movements directed toward warmer 
                           # temperature then current temperature 
                           # FALSE = random movements
                           mean.starting.depth = 20,  # starting depth 
                           sd.depth.excursion = 5,    # depth change allowance 
                           start.date = 273,          # Day of year of sart of simulation
                           end.date =  295,            # Day of year of end of simulation
                           deg.day.threshold = 0,# temperature threshold taken from Campbell, 1986
                           # used to calculate Perkins Eye Index
                           # or (egg developement) and
                           # Growing Degree Days
                           vertical.temperature.window = 220, # range over which temperature sampled 
                           vertical.movement.window = 15,   #max move range
                           temperature.selectivity.coefficient =  0, 
                           environment.data =  'Final .BayOfFundy.FVCOM.SeasonalBottomTemperature.RData'
) {
  n.days <- end.date - start.date   # number of days in simulation
  
  # read in the environmental data set here from FVCOM-GOM (NECOFS, CHEN et al, 2006)
  # load('Central.BayOfFundy.FVCOM.SeasonalBottomTemperature_1978_2014_Climatology.RData')
  #all figures done use data below
  #load('Final .BayOfFundy.FVCOM.SeasonalBottomTemperature.RData')
  load(environment.data)
  # create 1-D environment of mean temperature at depth for each day of the year
  depth <- depth.bins 
  temperature <- mean.seasonal.temperature.by.depth
  
  day.of.year <- seq(0,365,1)
  
  # allocate lobster state vectors
  lobster <- list(age = rep(0,n.lobsters),
                  lat = rep(0,n.lobsters),
                  lon = rep(0,n.lobsters),
                  depth = rep(0,n.lobsters),
                  current.temp = rep(0,n.lobsters),
                  max.temperature.target.depth = rep(0,n.lobsters), # max temp within informational range
                  max.temperature.target = rep(0,n.lobsters), # max temp within informational range
                  day.of.year = rep(0,n.lobsters),
                  PEI = rep(0,n.lobsters),            # Perkins Eye Index
                  pei.growth.rate = rep(0,n.lobsters),
                  cum.deg.days = rep(0,n.lobsters)    # cumulative degree days
  )
  
  # allocate lobster state time series
  lobster.archive <- list(age = array(0,dim=c(n.lobsters,n.days)),
                          lat = array(0,dim=c(n.lobsters,n.days)),
                          lon = array(0,dim=c(n.lobsters,n.days)),
                          depth = array(0,dim=c(n.lobsters,n.days)),
                          current.temp = array(0,dim=c(n.lobsters,n.days)),
                          max.temperature.target.depth = array(0,dim=c(n.lobsters,n.days)),
                          max.temperature.target = array(0,dim=c(n.lobsters,n.days)),
                          day.of.year = array(0,dim=c(n.lobsters,n.days)),
                          PEI = array(0,dim=c(n.lobsters,n.days)),
                          pei.growth.rate = array(0,dim=c(n.lobsters,n.days)),
                          cum.deg.days = array(0,dim=c(n.lobsters,n.days))
  )
  
  # start simulation
  current.day <- start.date
  
  # initialize the lobsters
  lobster$depth <- sample(mean.starting.depth+seq(-15,15,1),n.lobsters,replace=TRUE) 
  lobster$day.of.year <- rep(start.date,n.lobsters)
  lobster$current.temp <- temperature[cbind(match(lobster$depth,depth),current.day)]
  lobster$PEI<- 0 # would this be 0 because it is the initial state of the eggs?  
  lobster$pei.growth.rate<-0 
  lobster$cum.deg.days<- 0 #same here? 0 cum.deg.days ?
  
  # initialize data archive
  lobster.archive$depth[,1] <- lobster$depth
  lobster.archive$day.of.year[,1] <- lobster$day.of.year
  lobster.archive$current.temp[,1] <- lobster$current.temp
  lobster.archive$max.temperature.target.depth[,1] <- NA
  lobster.archive$max.temperature.target[,1] <- NA
  lobster.archive$PEI [,1]<-lobster$PEI
  lobster.archive$pei.growth.rate[,1]<-lobster$pei.growth.rate
  lobster.archive$cum.deg.days[,1]<- lobster$cum.deg.days
  
  
  for (d in seq(2,n.days,1)) {
    
    current.day <- current.day + 1
    # BEGIN update the state vector
    
    lobster$age <- lobster$age + 1  # track age
    lobster$day.of.year <- (lobster$day.of.year)%%365 + 1    # track day.of.year # modulo operator; after day 360 it starts at 1 again.    #lobster$day.of.year <- (lobster$day.of.year+1)%%365 + 1 - removed the +1 in bracket to make dayofyear one at a tim June 23 2020 ph
    lobster$current.temp <- pmax(0,temperature[cbind(match(lobster$depth,depth),current.day)])
    
    # temperature rule skews likelihood of lobster movement strongly to 
    # temperature warmer than its current temperature
    if (useTemperatureRule == TRUE){
      for (n in 1:n.lobsters) {
        
        # characterize daily environmental temperature profile relative to lobster's 
        # current temperature
        
        delta.T <-   temperature[,current.day] - lobster$current.temp[n]
        # global.delta.T <- delta.T
        
        # identify depth window of temperatures from which a lobster can choose
        # and ensure that scope remains within depth range of model domain
        vertical.temperature.scope <- seq(lobster$depth[n]-vertical.temperature.window,
                                          lobster$depth[n]+vertical.temperature.window,1)
        vertical.temperature.scope <- vertical.temperature.scope[which(vertical.temperature.scope > 0 & 
                                                                         vertical.temperature.scope < 220)]
        
        # potential range of movement
        # assumes that range of information for temperature equals or exceeds movement range
        # note that Pr(move) set to zero outside of movement range below
        vertical.scope <- vertical.temperature.scope
        
        # mask <- 0*depth
        # mask[vertical.temperature.scope] <- 1
        # identify depth window over which a lobster can move
        # and ensure that scope remains within depth range of model domain
        vertical.move.scope <- seq(lobster$depth[n]-vertical.movement.window,
                                   lobster$depth[n]+vertical.movement.window,1)
        vertical.move.scope <- vertical.move.scope[which(vertical.move.scope > 0 & 
                                                           vertical.move.scope < 250)]
        
        # extract relative temperature over informational depth range,scale to 0-1
        local.delta.T <- delta.T[vertical.temperature.scope]
        local.delta.T <- (local.delta.T-min(local.delta.T))/diff(range(local.delta.T))
        #   plot(local.delta.T,vertical.temperature.scope,ylim=c(max(vertical.temperature.scope),min(vertical.temperature.scope)), type = "l",lty = 1)
        #  plot(temperature[vertical.temperature.scope,current.day],vertical.temperature.scope,ylim=c(max(vertical.temperature.scope),min(vertical.temperature.scope)), xlim = c(4,14), type = "l", lwd = 2, ylab = "Depth (m)",xlab = "Bottom Temperature (C°)", las = 1)
        max.temperature.target.depth <- mean(vertical.temperature.scope[which(local.delta.T==max(local.delta.T))])
        lobster$max.temperature.target[n] <- max(local.delta.T)
        lobster$max.temperature.target.depth[n] <- max.temperature.target.depth
        
        
        # calculate localized, temperature-responsive movement       
        # uniform probability within specified depth range around current location
        base.move.llhd <- rep(0,length(vertical.temperature.scope))  
        base.move.llhd[vertical.temperature.scope%in%vertical.move.scope] <- 1
        base.move.llhd <- base.move.llhd/sum(base.move.llhd)
        
        # calculate probability of movement in response to temperature
        # exponential distribution, with greatest probability assigned to 
        # maximum temperature within "sensed" range; declining probability 
        # for increasingly cooler temperatures       
        # temperature.selectivity.coefficient <- 0.1  ### turn this into an argument; higher values increase selectivity for movement to higher temperatures
        # delta.depth.to.max.temp <- abs(max.temperature.target.depth - vertical.move.scope)
        delta.depth.to.max.temp <- abs(max.temperature.target.depth - vertical.temperature.scope)
        temperature.dependent.movement.tendency <- 
          exp(-delta.depth.to.max.temp * temperature.selectivity.coefficient)
        #plot(exp(-delta.depth.to.max.temp * temperature.selectivity.coefficient), type = "l", ylab = "Max Temperature Identified", xlab = "Depth (m)")
        #plot(delta.depth.to.max.temp, type = "l")
        #only keep values within movement range
        tdmt.tmp <- 0*temperature.dependent.movement.tendency
        tdmt.tmp[vertical.temperature.scope%in%vertical.move.scope] <- 
          temperature.dependent.movement.tendency[vertical.temperature.scope%in%vertical.move.scope]
        temperature.dependent.movement.tendency <- tdmt.tmp
        
        temperature.dependent.movement.tendency <- temperature.dependent.movement.tendency/sum(temperature.dependent.movement.tendency)
        #plot(temperature.dependent.movement.tendency, type ="l", ylab = "Movement Tendency", xlab = "Depth (m)", cex.lab = 1.25, cex.axis = 1.25, ylim = c(0,0.6) )
        delta.T.quasi.llhd <- temperature.dependent.movement.tendency
        
        # delta.T.quasi.llhd <- exp(-temperature.selectivity.coefficient*(1 - local.delta.T))  # coefficient affects concentration of probability towards max(T) (higher = more concentrated/steep)
        # delta.T.quasi.llhd <- delta.T.quasi.llhd/sum(delta.T.quasi.llhd)
        # need to map temperature profile to movement range.
        
        # combine depth-dependent and temperature-dependent Pr(movement)
        tdep.move.llhd <- base.move.llhd * temperature.dependent.movement.tendency
        #plot(base.move.llhd)
        tdep.move.llhd <- tdep.move.llhd/sum(tdep.move.llhd)
        # plot(tdep.move.llhd,vertical.temperature.scope,ylim=c(max(vertical.temperature.scope),min(vertical.temperature.scope)), type= "l") #probability of moving to depth
        #         points(base.move.llhd,vertical.temperature.scope,col='red')
        move.llhd <- tdep.move.llhd
        move.llhd.vector <- cumsum(move.llhd/sum(move.llhd)) # convert to cumulative probability distribution
        
        # randomly select position along cumulative probability vector
        draw <- runif(1,0,1)  ### further opportunity here to vectorize?!? 
        move.idx <- min(which(move.llhd.vector >= draw),na.rm=TRUE)
        
        lobster$depth[n] <- pmax(10,pmin(250,vertical.scope[move.idx]))
        # print(lobster$depth[n])
      }
      # else the lobster moves randomly up or down 5 meters depth from current location
    } else {
      lobster$depth <- pmax(10,  # elementwise comparison (max)
                            pmin(220,
                                 round(lobster$depth + runif(n = n.lobsters,min = -15,max = 15))))
      # print(lobster$depth)
    }
    # pmax (0) to ensure no negative temperatures which blow up pei.growth.rate formula
    # tracking embryo development in relation to temperature exposure
    # associated with lobster trajectory
    lobster$pei.growth.rate <- 1000*(0.00002*lobster$current.temp^2.25008) # 1000 is to convert mm to um
    lobster$PEI <- lobster$PEI + lobster$pei.growth.rate
    lobster$cum.deg.days <- lobster$cum.deg.days + pmax(0,lobster$current.temp - deg.day.threshold)  
    
    
    # update the archive 
    lobster.archive$age[,d] <- lobster$age
    lobster.archive$day.of.year[,d] <- lobster$day.of.year
    lobster.archive$depth[,d] <- lobster$depth
    lobster.archive$current.temp[,d] <- lobster$current.temp
    lobster.archive$max.temperature.target.depth[,d] <- lobster$max.temperature.target.depth
    lobster.archive$max.temperature.target[,d] <- lobster$max.temperature.target
    lobster.archive$PEI[,d]<-lobster$PEI
    lobster.archive$pei.growth.rate[,d]<-lobster$pei.growth.rate
    lobster.archive$cum.deg.days[,d]<-lobster$cum.deg.days
  }
  
  
  Lobster.Simulation <- list(lobster.archive = lobster.archive,
                             start.date = start.date,
                             end.date = end.date,
                             useTemperatureRule = useTemperatureRule)
  
  return(Lobster.Simulation)
  
  
  
}



####  ARCHIVED alternate formulations of movement rules
# delta.T.quasi.llhd <- exp(-delta.T)
# delta.T.quasi.llhd <- delta.T.quasi.llhd/sum(delta.T.quasi.llhd)
# delta.T.quasi.llhd[delta.T.quasi.llhd < 1] <- delta.T.quasi.llhd[delta.T.quasi.llhd < 1]/10
# move.llhd <- dnorm(depth,lobster$depth[n],sd.depth.excursion)
# move.llhd <- dunif(depth,lobster$depth[n]-5,lobster$depth[n]+5)
# move.llhd <-  move.llhd * delta.T.quasi.llhd
# move.llhd <-  delta.T.quasi.llhd


# base.move.llhd <- dnorm(vertical.temperature.scope,lobster$depth[n],sd.depth.excursion) # random Normal variate
