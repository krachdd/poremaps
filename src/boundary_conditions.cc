/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright 2024 David Krach, Matthias Ruf

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the “Software”), to deal in 
the Software without restriction, including without limitation the rights to use, 
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.


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
    int j, k;
    int zmax_procs = dims[2] - 1; // index of last rank in the z direction
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