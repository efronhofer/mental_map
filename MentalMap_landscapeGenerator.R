##########################################################################
# Mental map
#
# Emanuel A. Fronhofer et al.
#
# September 2013 -- landscape generator
###########################################################################

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
# FUNCTIONS

# generate key to compare two matrices
key <- function(m){
  paste(sep="\1",m[,1],m[,2])
}

# simulate hexagonal grid for plotting results
hexsim <- function(xo,yo){
  if (yo %% 2 == 0) {
		xn <- xo + 0.5
	}else{
		xn <- xo + 0.0
	}
	return(xn)
}

##########################################################################
# INITIALIZATION

# location of all patches in hexagnal grid
world[,1] <- rep(c(1:DIM,1:(DIM-1)),DIM/2)
world[,2] <- rep(c(1:DIM),rep(c(DIM,DIM-1),DIM/2))
world[,1] <- mapply(hexsim,world[,1],world[,2])
# correct exact y position in the hexagonal grid
world[,2] <- world[,2]* sqrt(3)/2

world[,3] <- 0

###########################################################################
# add clumped resources (Thomas process)
# note: rThomas(number of clusters, displacement from cluster center, points per cluster)
# less dense patches
#patches8 <- unique(round(as.data.frame(rThomas(DIM/2,3/DIM,DIM))*DIM),digits=1)
#dense patches
# dense
if (landscape=="dense"){patches8 <- unique(round(as.data.frame(rThomas(0.005*DIM^2,1/DIM,150))*DIM))}
#if (landscape=="dense"){patches8 <- unique(round(as.data.frame(rThomas(DIM^2/200,0.0075,200))*DIM))}

# intermediate
if (landscape=="intermediate"){ patches8 <- unique(round(as.data.frame(rThomas(0.005*DIM^2,2/DIM,35))*DIM))}
#if (landscape=="intermediate"){ patches8 <- unique(round(as.data.frame(rThomas(DIM^2/200,0.025,50))*DIM))}
  
# random
if (landscape=="random"){ patches8 <- unique(round(as.data.frame(rThomas(0.005*DIM^2,10/DIM,26))*DIM))}
 
patches <- cbind(mapply(hexsim,patches8[,1],patches8[,2]),patches8[,2]*sqrt(3)/2)

# save resources to world
world[which(is.element(key(world[,1:2]),key(patches))),3] <- 1
world[,4] <- world[,3]

# calculate E0 (mean probability of resource encounter)
E0 <- length(patches)/length(world[,1])

##################################################################################
