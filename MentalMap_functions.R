##########################################################################
# Mental map
#
# Emanuel A. Fronhofer et al.
#
# September 2013 -- functions
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
# FUNCTIONS

# movement in a hexagonal grid
newloc <- function(pos, nd){
  if (nd == 1) npos <- pos + (DIM-1)
	if (nd == 2) npos <- pos + DIM
	if (nd == 3) npos <- pos + 1
	if (nd == 4) npos <- pos - (DIM-1)
	if (nd == 5) npos <- pos - DIM
	if (nd == 6) npos <- pos - 1
	return(npos)
}

findLine <- function(noR,nr,pos){

      a1 <- (pos-nr+noR*(-DIM+1))
      a2 <- (pos-nr+noR*DIM)
      b2 <- pos+nr+noR*(DIM-1)
      b1 <- pos+nr-noR*DIM
      
      return(c(c(a1:b1), c(a2:b2)))
}

    
findRings <- function(pos,nr){
  if(nr==1){
    rc <- sapply(1:6,newloc,pos=pos)
  }else{

    rc <- c( c((pos-nr):(pos+nr))[-(nr+1)] , unlist(sapply(1:nr,findLine,nr=nr,pos=pos)) )
        
  }  
  return(rc)
}


# find nearest cell which contained a resource item am return distance
resProb <- function(actcell1, II){
  
  # check whether something is known about this cell
  E_prel <- map[[II]][actcell1,3]
  
  if(is.na(E_prel)){
    if (nRings > 0){
    # find all cells in nRings rings
    #nnc <- sapply(1:6,newloc,pos=actcell1)
    nnc <- findRings(actcell1,nRings)
    # save info
    info <- map[[indiv]][nnc,4]
    info[which(is.na(info))] <- E0
    
    a <- mean(info)
    }else{
      a <- E0
    }
    
  }else{
    a <- E_prel
  }
  
  return(a)
  
}


# calculate attractivity of taking a certain direction from mental map
# including only one additional ring of cells
E_Pmap <- function(actcell, d, II){
  # determine new cell from act cell and direction
	nc <- newloc(actcell,d)
	# save attractivity from map
	E_newcell <- map[[II]][nc,3]
	# determine all cell
	dirs <- c(d-1,d,d+1)
	dirs[which(dirs == 0)] <- 6
	dirs[which(dirs == 7)] <- 1
	sc <- sapply(dirs,newloc,pos=nc)
	# save attractivities from map
  E_neighbours <- sapply(sc,resProb,II=II)
	#E_neighbours <- sapply(sc,resProb, II=II)
  nextcells <- c( newloc(sc[1],dirs[1]), sapply(dirs,newloc,pos=sc[2])  , newloc(sc[3],dirs[3]) )
  E_nextNeighbours <- sapply(nextcells,resProb,II=II)

	attract <- sum(E_newcell, E_neighbours, E_nextNeighbours)
  #attract <- sum(E_newcell, E_neighbours)

	return(attract)

}


# calculate attractivity of taking a certain direction from mental map
# including only one additional ring of cells
# only the two rings in the perc range
E_Pmap1 <- function(actcell, d, II){
  # determine new cell from act cell and direction
	nc <- newloc(actcell,d)
	# save attractivity from map
	E_newcell <- map[[II]][nc,3]
	# determine all cell
	dirs <- c(d-1,d,d+1)
	dirs[which(dirs == 0)] <- 6
	dirs[which(dirs == 7)] <- 1
	sc <- sapply(dirs,newloc,pos=nc)
	# save attractivities from map
	E_neighbours <- sapply(sc,resProb, II=II)

	attract <- sum(E_newcell, E_neighbours)

	return(attract)

}
