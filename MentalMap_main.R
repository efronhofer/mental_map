##########################################################################
# Mental map
#
# Emanuel A. Fronhofer
#
# September 2013
##########################################################################

# Copyright (C) 2013  Emanuel A. Fronhofer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

##########################################################################

# Instead of perceiving/ knowing the location of every resource patch in
# the world like in the biased correlated random walk, we assume that indi-
# viduals gain information during displacement. This information is saved
# as a path vector. From this information the animal constructs a
# probability map of findig resources.
# If no resources are found the animal tries to maximise the amount of
# information it gains since the mean probability of findig a resource 
# is >0 in unkown cells.
#
# In this version absorbing boundary conditions are assumend, memory is 
# limited in time (forgetting!), the animal is spatially informed, i.e.
# it incorporates assumptions about the envirnoment, i.e. spatial 
# correlation of resources into its cognitive map.
#
# This version contains a possibility to introduce "reading errors" while 
# retreiving information from the mental map (parameter epsilon), more
# than one individual and regrowing resources (note that this only works 
# with one individual).
#
# if more than one direction is equally attractive up to 3 rings of the 
# mental map are used

#########################################################################
# clear workspace
rm(list=ls())

# libraries for generating landscape (Thomas process)
library(mgcv)
library(spatstat)

# library for plotting of mental map
library(akima)

#########################################################################
# PARAMETERS
source("MentalMap_parameters.R")

########################################################################
# FUNCTIONS
source("MentalMap_functions.R")

##########################################################################
# DATA STRUCTURES

# saves x and y coordinates and resource content
world <- matrix(ncol=4,nrow=DIM^2-DIM/2)

# saves x and y coordinates; actual resources(after depeltion), observed resoures, time of observation
map <- list()
for (indiv in 1:N_Ind){map[[indiv]] <- matrix(ncol=5,nrow=DIM^2-DIM/2)}

# saves patch address
path <- list()
for (indiv in 1:N_Ind){path[[indiv]] <- numeric(TMAX)}

# saves patch location
path_coords <- list()
for (indiv in 1:N_Ind){path_coords[[indiv]] <- matrix(ncol=2,nrow=TMAX)}

#########################################################################
# INITIALIZATION

# landscape generation
source("MentalMap_landscapeGenerator.R")

# intialize map
for (indiv in 1:N_Ind){
  map[[indiv]][,1:2] <- world[,1:2]
  map[[indiv]][,3] <- NA
  map[[indiv]][,4] <- NA
  map[[indiv]][,5] <- NA
}

# cex for plotting
if (DIM == 100) {
  dotsize <- 0.8
  }else{
  dotsize <- 0.35
  }

# plot world
x11(width=12,height=6)
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
plot(seq(1,DIM,len=100),seq(1,DIM*sqrt(3)/2,len=100),type="n",xaxt="n",yaxt="n",xlab="n",ylab="n")
points(world[which(world[,3]==(1)),1],world[which(world[,3]==(1)),2], col="green3",pch=16, cex=dotsize)

# starting patch (center of world)
center <- round(DIM*DIM/2- DIM/4 - DIM/2)

# start in the centre if there is only one individual
if (N_Ind == 1){
  path[[indiv]][1] <- center
  path_coords[[indiv]][1,1] <- world[path[[indiv]][1],1]
  path_coords[[indiv]][1,2] <- world[path[[indiv]][1],2]
}else{
  # randomly assign starting point if there are more individuals
  for (indiv in 1:N_Ind){
    path[[indiv]][1] <- center - sample(c((-DIM/4):(DIM/4)),1) - DIM*sample(c((-DIM/4):(DIM/4)),1)
    path_coords[[indiv]][1,1] <- world[path[[indiv]][1],1]
    path_coords[[indiv]][1,2] <- world[path[[indiv]][1],2]
  }
}

# counter for no. of patches found
pfound <- rep(0,N_Ind)

# save indiv timestep
timestep <- numeric(N_Ind)

# save list of depleted resources
depResList <- numeric(0)

# calculate neg. exp. forgetting
#forgetFact <- exp(-f)

# init direction
olddir <- sample(1:6,1)

#########################################################################
# ITERATION

for (iteration in 1:(TMAX-1)){
    
    # once per timestep: forgetting
    forgottenCells <- which(map[[indiv]][,5] < (iteration-forgetLag))
    map[[indiv]][forgottenCells,3:5] <- NA
    
    # randomly choose individual
    indiv <- sample(1:N_Ind,1)
  
    # indiv spec timestep (counter for path)
    timestep[indiv] <- timestep[indiv] +1
  
  	# plot path
    lines(path_coords[[indiv]][,1],path_coords[[indiv]][,2])
    
    # resource regrowth
    if (length(depResList) >= growthLag){
       # create new resource item
      world[depResList[iteration-growthLag],3] <- E0 + 1
      # plot
      points(world[depResList[iteration-growthLag],1],world[depResList[iteration-growthLag],2],pch=16,col="green3",cex=dotsize)
    }
  
  	# save actual amount of resources found
  	ares <- world[path[[indiv]][timestep[indiv]],3]
  
  	# check whether there was a patch, deplete it and plot
  	if (ares != 0){
  		# deplete resources
  		world[path[[indiv]][timestep[indiv]],3] <- 0
      # save patch no to list
      depResList[length(depResList)+1] <- path[[indiv]][timestep[indiv]]
  		# plot
  		points(path_coords[[indiv]][timestep[indiv],1], path_coords[[indiv]][timestep[indiv],2],pch=16,col="red",cex=dotsize)
  		# increase counter for no. of patches found
  		pfound[indiv] <- pfound[indiv] + 1
  	}else{
      depResList[length(depResList)+1] <- NA
  	}
  	
    # incorporate perceptual range
    # find the six neighbouring cells
    nn6 <- sapply(c(1:6),newloc,pos=path[[indiv]][timestep[indiv]])

    if (is.na(map[[indiv]][path[[indiv]][timestep[indiv]],3])){
    # hier muss rein: nur überschreiben, wenn unbekannt, v.a. [,4]!
    # save perceived information to map
  	# 1. actual cell
    map[[indiv]][path[[indiv]][timestep[indiv]],3] <- 0
    map[[indiv]][path[[indiv]][timestep[indiv]],4] <-  world[path[[indiv]][timestep[indiv]],3]
    map[[indiv]][path[[indiv]][timestep[indiv]],5] <- iteration
    }else{
      map[[indiv]][path[[indiv]][timestep[indiv]],3] <- world[path[[indiv]][timestep[indiv]],3]
      map[[indiv]][path[[indiv]][timestep[indiv]],5] <- iteration
    }
    
    # 2. six nearest neighbours (perceptual range)
    for (help in 1:6){
      if (is.na(map[[indiv]][nn6[help],3])){
      map[[indiv]][nn6[help],3] <- world[nn6[help],3]
      map[[indiv]][nn6[help],4] <- world[nn6[help],3]
      map[[indiv]][nn6[help],5] <- iteration
      }else{
        map[[indiv]][nn6[help],3] <- world[nn6[help],3]
        map[[indiv]][nn6[help],5] <- iteration
      }
    }

  	
    # check back direction
    backdir <- olddir+3
    if (backdir > 6) backdir <- backdir-6
    # possible directions for movement (all except back)
    posdirs <- c(1:6)[-backdir]    
    
    # check only nn6
    EW_nn6 <-  map[[indiv]][nn6[posdirs],3]
    
    direction <- posdirs[which(EW_nn6==max(EW_nn6))]
    
    if (length(direction) > 1){
      
      # look into map and calc. attractivity of all directions
      EW_1 <- sapply(direction,E_Pmap1,actcell=path[[indiv]][timestep[indiv]],II=indiv)
      
      # take most attractive route
      # if more than one option is possible, keep old direction or throw a dice!
      direction <- direction[which(EW_1==max(EW_1))]
      if (length(direction) > 1){
      
        # look into map and calc. attractivity of all directions
        EW_0 <- sapply(direction,E_Pmap,actcell=path[[indiv]][timestep[indiv]],II=indiv)
       
        # introduce some stochasticity during retreiving information
        EW <- EW_0 + rnorm(length(EW_0),0,epsilon)
        # round in order to dismiss error of hex grid (sqrt(3)/2)
        EW <- round(EW,digits=10)
       
        # take most attractive route
        # if more than one option is possible, keep old direction or throw a dice!
      	direction <- direction[which(EW==max(EW))]
      	if (length(direction) > 1){
          if (is.element(olddir, direction)){
              if (fsl == 1){
                direction <- olddir
              }else{
                direction <- sample(direction,1)
              }
            }else{
              direction <- sample(direction,1)
            }
      	}
    
    }
    }
    olddir <- direction

    #print(direction)
    
  	path[[indiv]][timestep[indiv]+1] <- newloc(path[[indiv]][timestep[indiv]],direction)
  	path_coords[[indiv]][timestep[indiv]+1,1] <- world[path[[indiv]][timestep[indiv]+1],1]
  	path_coords[[indiv]][timestep[indiv]+1,2] <- world[path[[indiv]][timestep[indiv]+1],2]
  
    # check for x-cood
    if (is.element(round(path_coords[[indiv]][timestep[indiv]+1,1],digits=4),c(1,1.5,DIM-1,DIM-0.5,DIM))) break   
    # check for y-coord	
    if (is.element(round(path_coords[[indiv]][timestep[indiv]+1,2],digits=4),round(c(1:(nRings+3),DIM:(DIM-(nRings+3)))*sqrt(3)/2,digits=4))) break
    
}

########################################################################
# return success of searching strategy
print(paste('no. of patches found: ',pfound))
print(paste('no. of timesteps: ',timestep))
print(paste('detection rate: ',pfound/timestep))

########################################################################

# plot mental map
plotmap <- map[[1]][,3] + map[[1]][,4]

#zz <- interp(map[[1]][,1],map[[1]][,2],plotmap,seq(1,DIM,length=DIM),seq(1,DIM,length=DIM))

par(mar=c(1,1,1,1))
plot(seq(1,DIM,len=100),seq(1,DIM*sqrt(3)/2,len=100),type="n",xaxt="n",yaxt="n",xlab="n",ylab="n")
#image(zz,zlim=c(0,2),add=T,col=c("white","grey","green","red"))
points(map[[1]][,2][which(plotmap==2)]~map[[1]][,1][which(plotmap==2)],pch=16, cex=dotsize, col="green")
points(map[[1]][,2][which(plotmap==1)]~map[[1]][,1][which(plotmap==1)],pch=16, cex=dotsize, col="red")
points(map[[1]][,2][which(plotmap==0)]~map[[1]][,1][which(plotmap==0)],pch=16, cex=dotsize, col="grey")

#grey(seq(1,0,len=4))

# save data for MS plot
write.table(patches, file="./Data/AllPatches_T0.out",row.names=F,col.names=F)
write.table(depResList, file="./Data/DepletedResources.out",row.names=F,col.names=F)
write.table(map, file="./Data/Map_Tend.out",row.names=F,col.names=F)
write.table(path_coords[[indiv]], file="./Data/Path.out",row.names=F,col.names=F)
write.table(world, file="./Data/World_Tend.out",row.names=F,col.names=F)


