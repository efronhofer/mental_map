##########################################################################
# Mental map
#
# Emanuel A. Fronhofer et al.
#
# September 2013 -- parameter file
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
# PARAMETERS

# world dimensions
DIM <- 100

# simulation time (number of steps)
TMAX <- 1000

# number of animals
N_Ind <- 1

# no of rings
nRings <- 2

# stochasticity when retreiving infromation (sd of gaussion dist around 1)
epsilon <- 0

# landscape configuration, i.e. clumping (one of "dense", "intermediate", "random")
landscape <- "dense"

# lag (timesteps) for resource regrowth
growthLag <- TMAX

# force straight line
fsl <- 0

# foretting time lag
forgetLag <- TMAX
