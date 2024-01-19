/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright (C) 2023  David Krach, Matthias Ruf

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


 * Authors:    David Krach 
 * Date:       2021/12
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    void reapply_pressure_gradient:
 *             Re-Set pressure Dirichlet BC
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/

// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "constants.h"
// **********************************************************************

void reapply_pressure_gradient(bool*** proc_geom, 
                               int* size, 
                               int* new_proc_size, 
                               double*** press,
                               MPI_Comm comm_cart, 
                               int* cart_coords, 
                               int* dims)
{
    int i, j, k;
    int zmax_procs = dims[2] - 1; // max number of process in z dir in cart communicator
    int lim = 2;
    if (cart_coords[2] == 0){ // reset inlet pressure 
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                if (proc_geom[1][j][k] == 0){
                    press[1][j][k] = size[2] + 2; // left boundary
                }
            }
        }
    }

    if (cart_coords[2] == zmax_procs){ // reset outlet pressure
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                if (proc_geom[new_proc_size[2]-2][j][k] == 0){
                    press[new_proc_size[2]-2][j][k] = 2.0; // right boundary
                }
                if (proc_geom[new_proc_size[2]-3][j][k] == 0){
                    press[new_proc_size[2]-3][j][k] = 3.0; // right boundary
                }
            }
        }
    }
}