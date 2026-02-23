/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright 2024 David Krach, Matthias Ruf

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
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
 * Purpose:    void allocate_initial_field:
 *             Allocate and initialise velocity and pressure fields,
 *             then communicate pressure halos.
 *
 *             void determine_comm_pattern:
 *             Set per-rank bc[] array that controls which directions
 *             communicate (MPI) and which apply a wall condition.
 *
 *             void set_halos_initial:
 *             Zero velocity halo layers for no-slip boundary ranks.
 * Contents:
 *
 * Changelog:  Version: 1.0.0
 *
 * Each rank stores bc[3], one entry per spatial direction [x, y, z].
 * Values encode which wall (if any) the rank borders and the wall type:
 *
 *   0:  communicate in both directions (interior rank)
 *   1:  communicate in +1, slip condition in -1
 *   2:  communicate in +1, no-slip condition in -1
 *   3:  communicate in -1, slip condition in +1
 *   4:  communicate in -1, no-slip condition in +1
 *   5:  slip condition in both -1 and +1
 *   6:  no-slip condition in both -1 and +1
 *
 ************************************************************************/


// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include "mpi.h"
#include "constants.h"
#include "parallelization.h"
// **********************************************************************


void allocate_initial_field(bool*** proc_geom,
                            int* size,
                            int* new_proc_size,
                            double*** press,
                            double*** vel_x,
                            double*** vel_y,
                            double*** vel_z,
                            MPI_Comm comm_cart,
                            int* neighbors,
                            int* starts,
                            unsigned int* bc)
{
    int i, j, k;

    // Allocate velocity and pressure arrays; set all fields to zero/initial values.
    for (i = 0; i < new_proc_size[2]; i++){
        vel_x[i] = (double**)malloc(new_proc_size[1]*sizeof(double*));
        vel_y[i] = (double**)malloc(new_proc_size[1]*sizeof(double*));
        vel_z[i] = (double**)malloc(new_proc_size[1]*sizeof(double*));
        press[i] = (double**)malloc(new_proc_size[1]*sizeof(double*));
        for (j = 0; j < new_proc_size[1]; j++){
            vel_x[i][j] = (double*)malloc(new_proc_size[0]*sizeof(double));
            vel_y[i][j] = (double*)malloc(new_proc_size[0]*sizeof(double));
            vel_z[i][j] = (double*)malloc(new_proc_size[0]*sizeof(double));
            press[i][j] = (double*)malloc(new_proc_size[0]*sizeof(double));
            for (k = 0; k < new_proc_size[0]; k++){
                vel_x[i][j][k] = 0.0;
                vel_y[i][j][k] = 0.0;
                vel_z[i][j][k] = 0.0;
                press[i][j][k] = -1.0;
            }
        }
    }

    // Initialise pressure with a linear gradient along z (flow direction).
    // Fluid voxels get a value proportional to their z position; solid voxels get 0.
    int lim = 2;
    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                if (proc_geom[i][j][k] == 0){
                    press[i][j][k] = (size[2] - starts[2] - i) + lim + 2;
                }
                else {
                    press[i][j][k] = 0.0;
                }
            }
        }
    }

    // Communicate pressure values into halo layers.
    communicate_halos(press, new_proc_size, neighbors, MPI_DOUBLE, comm_cart, bc, 0);
}


void determine_comm_pattern(int* dims,
                            int rank,
                            int* cart_coords,
                            unsigned int* bc,
                            int bc_method,
                            int* neighbors)
{
/*
 * bc_method selects the physical setup:
 *
 *   0:  fully periodic in all directions
 *       → all ranks communicate in all directions
 *
 *   1:  periodic in z; slip walls in x and y
 *       → all ranks communicate in z; boundary ranks in x/y apply slip
 *
 *   2:  periodic in z; no-slip walls in x and y
 *       → all ranks communicate in z; boundary ranks in x/y apply no-slip
 *
 *   3:  non-periodic; slip walls in x and y; pressure-driven flow in z
 *       → boundary ranks in all directions apply slip
 *
 *   4:  non-periodic; no-slip walls in x and y; pressure-driven flow in z
 *       → boundary ranks in x/y apply no-slip; boundary ranks in z apply slip
 */

    int i;

    // Default: all ranks communicate in both directions (interior behaviour).
    // Boundary ranks below override their relevant entry.
    for (i = 0; i < ndims; i++){ bc[i] = 0; }

    if (bc_method == 0){
        // Fully periodic — bc stays 0 for all ranks and all directions.
    }

    else if (bc_method == 1){
        // Periodic z; slip x and y walls.
        bc[2] = 0;  // z is periodic: all ranks communicate

        if (dims[0] == 1){ bc[0] = 5; }          // single rank in x: slip on both sides
        else if (dims[0] > 1){
            if (cart_coords[0] == 0)           { bc[0] = 1; }  // first rank in x: slip at -x
            if (cart_coords[0] + 1 == dims[0]) { bc[0] = 3; }  // last rank in x:  slip at +x
        }

        if (dims[1] == 1){ bc[1] = 5; }
        else if (dims[1] > 1){
            if (cart_coords[1] == 0)           { bc[1] = 1; }
            if (cart_coords[1] + 1 == dims[1]) { bc[1] = 3; }
        }
    }

    else if (bc_method == 2){
        // Periodic z; no-slip x and y walls.
        bc[2] = 0;

        if (dims[0] == 1){ bc[0] = 6; }
        else if (dims[0] > 1){
            if (cart_coords[0] == 0)           { bc[0] = 2; }
            if (cart_coords[0] + 1 == dims[0]) { bc[0] = 4; }
        }

        if (dims[1] == 1){ bc[1] = 6; }
        else if (dims[1] > 1){
            if (cart_coords[1] == 0)           { bc[1] = 2; }
            if (cart_coords[1] + 1 == dims[1]) { bc[1] = 4; }
        }
    }

    else if (bc_method == 3){
        // Non-periodic; slip walls in all directions.
        if (dims[0] == 1){ bc[0] = 5; }
        else if (dims[0] > 1){
            if (cart_coords[0] == 0)           { bc[0] = 1; }
            if (cart_coords[0] + 1 == dims[0]) { bc[0] = 3; }
        }

        if (dims[1] == 1){ bc[1] = 5; }
        else if (dims[1] > 1){
            if (cart_coords[1] == 0)           { bc[1] = 1; }
            if (cart_coords[1] + 1 == dims[1]) { bc[1] = 3; }
        }

        if (dims[2] == 1){ bc[2] = 5; }
        else if (dims[2] > 1){
            if (cart_coords[2] == 0)           { bc[2] = 1; }
            if (cart_coords[2] + 1 == dims[2]) { bc[2] = 3; }
        }
    }

    else if (bc_method == 4){
        // Non-periodic; no-slip x and y walls; slip z walls (pressure-driven).
        if (dims[0] == 1){ bc[0] = 6; }
        else if (dims[0] > 1){
            if (cart_coords[0] == 0)           { bc[0] = 2; }
            if (cart_coords[0] + 1 == dims[0]) { bc[0] = 4; }
        }

        if (dims[1] == 1){ bc[1] = 6; }
        else if (dims[1] > 1){
            if (cart_coords[1] == 0)           { bc[1] = 2; }
            if (cart_coords[1] + 1 == dims[1]) { bc[1] = 4; }
        }

        // z direction uses slip (not no-slip) even for bc_method 4 — the pressure
        // gradient is imposed via reapply_pressure_gradient, not a wall condition.
        if (dims[2] == 1){ bc[2] = 5; }
        else if (dims[2] > 1){
            if (cart_coords[2] == 0)           { bc[2] = 2; }
            if (cart_coords[2] + 1 == dims[2]) { bc[2] = 4; }
        }
    }
}


void set_halos_initial( double*** var,
                        int* new_proc_size,
                        unsigned int* bc)
{
    // Zero the halo layers that correspond to no-slip wall boundaries.
    // Only called for velocity fields; pressure halos are handled separately.
    int i, j;

    // X halos: indices 0,1 (left wall) and nx-2,nx-1 (right wall).
    for (i = 0; i < new_proc_size[2]; i++){
        for (j = 0; j < new_proc_size[1]; j++){
            if (bc[0] == 2 || bc[0] == 6){
                var[i][j][0] = 0.0;
                var[i][j][1] = 0.0;  // fixed: was erroneously set to [0] twice
            }
            if (bc[0] == 4 || bc[0] == 6){
                var[i][j][new_proc_size[0]-1] = 0.0;
                var[i][j][new_proc_size[0]-2] = 0.0;
            }
        }
    }

    // Y halos: indices 0,1 (bottom wall) and ny-2,ny-1 (top wall).
    for (i = 0; i < new_proc_size[2]; i++){
        for (j = 0; j < new_proc_size[0]; j++){
            if (bc[1] == 2 || bc[1] == 6){
                var[i][0][j] = 0.0;
                var[i][1][j] = 0.0;
            }
            if (bc[1] == 4 || bc[1] == 6){
                var[i][new_proc_size[1]-1][j] = 0.0;
                var[i][new_proc_size[1]-2][j] = 0.0;
            }
        }
    }

    // Z halos: indices 0,1 (front wall) and nz-2,nz-1 (back wall).
    for (i = 0; i < new_proc_size[1]; i++){
        for (j = 0; j < new_proc_size[0]; j++){
            if (bc[2] == 2 || bc[2] == 6){  // fixed: was bc[3] (out of bounds)
                var[0][i][j] = 0.0;
                var[1][i][j] = 0.0;
            }
            if (bc[2] == 4 || bc[2] == 6){  // fixed: was bc[3] (out of bounds)
                var[new_proc_size[2]-1][i][j] = 0.0;
                var[new_proc_size[2]-2][i][j] = 0.0;
            }
        }
    }
}
