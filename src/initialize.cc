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
 * Purpose:    void allocate_initial_field:
 *             Initialize fields and commuincate halos
 * 
 *             void determine_comm_pattern:
 *             Get ranks communication pattern depending on 
 *             domain decomposition and BCs
 *
 *             void set_halos_initial:
 *             Set vel values in halos
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0

Each rank stores a array with 3 ints with holds information on how to 
communicate in each spatial direction [x, y, z] = {0, 1, ..., 6}

options:    0:  communicate in both directions
            1:  communicate in +1, slip condition in -1 
            2:  communicate in +1, no-slip condition in -1
            3:  communicate in -1, slip condition in +1 
            4:  communicate in -1, no-slip condition in +1
            5:  slip condition in -1 and +1 
            6:  no-slip condition in -1 and +1 

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
    int    i, j , k; 
    double gradp;

    // allocate velocities and set initial value 
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
                if (proc_geom[i][j][k] == 0){ // initiate velocity for fluid voxels 
                    vel_z[i][j][k] = 0.0; 
                }
                else {
                    vel_z[i][j][k] = 0.0;
                }
                vel_x[i][j][k] = 0.0;
                vel_y[i][j][k] = 0.0;
                press[i][j][k] = -1.0;

            }
        }
    }
    // initialize pressure gradient
    int lim = 2;
    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                // distinguish between solid and fluid voxels
                if (proc_geom[i][j][k] == 0){
                    // press[i][j][k] = (starts[2]+i-lim)+2;
                    press[i][j][k] = (size[2]-starts[2]-i)+lim+2;
                }
                else {
                    press[i][j][k] = 0.0;
                }
            }
        }
    }
    // MPI_Datatype etype = MPI_DOUBLE;
    // Communicate pressure values in halos
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
Possible input options in bc_method

options:    0:  requires:       domain periodic in all spatial directions
                communication:  all ranks in all directions 
                set:            nothing

            1:  requires:       domain periodic in z direction
                communication:  all ranks in z directions , non-boundary interfaces 
                set:            slip condition on x,y-boundary

            2:  requires:       domain periodic in z direction
                communication:  all ranks in z directions , non-boundary interfaces 
                set:            no-slip condition on x,y-boundary 

            3:  requires:       nothing
                communication:  non-boundary interfaces 
                set:            slip condition on x,y-boundary , continouus flow in z direction

            4:  requires:       nothing
                communication:  non-boundary interfaces 
                set:            no-slip condition on x,y-boundary , continouus flow in z direction

These are the 7 options to determine dependend on case and position of the rank in cart_coord

options:    0:  communicate in both directions
            1:  communicate in +1, slip condition in -1 
            2:  communicate in +1, no-slip condition in -1
            3:  communicate in -1, slip condition in +1 
            4:  communicate in -1, no-slip condition in +1
            5:  slip condition in -1 and +1 
            6:  no-slip condition in -1 and +1 
*/

    int i, j, k;

    if (bc_method == 0){
        for (i = 0; i < ndims; i++){bc[i] = 0;}
    }
    
    else if (bc_method == 1){
        bc[2] = 0;
        if (dims[0] == 1){ bc[0] = 5;}
        else if (dims[0] > 1){
            if (cart_coords[0] == 0){            bc[0] = 1;}
            if (cart_coords[0] + 1 == dims[0]){  bc[0] = 3;}
        }

        if (dims[1] == 1){ bc[1] = 5;}
        else if (dims[1] > 1){
            if (cart_coords[1] == 0){            bc[1] = 1;}
            if (cart_coords[1] + 1 == dims[1]){  bc[1] = 3;}

        }
    }

    else if (bc_method == 2){
        bc[2] = 0;
        if (dims[0] == 1){ bc[0] = 6;}
        else if (dims[0] > 1){
            if (cart_coords[0] == 0){            bc[0] = 2;}
            if (cart_coords[0] + 1 == dims[0]){  bc[0] = 4;}
        }

        if (dims[1] == 1){ bc[1] = 6;}
        else if (dims[1] > 1){
            if (cart_coords[1] == 0){            bc[1] = 2;}
            if (cart_coords[1] + 1 == dims[1]){  bc[1] = 4;}

        }
    }

    else if (bc_method == 3){
        if (dims[0] == 1){ bc[0] = 5;}
        else if (dims[0] > 1){
            if (cart_coords[0] == 0){            bc[0] = 1;}
            if (cart_coords[0] + 1 == dims[0]){  bc[0] = 3;}
        }

        if (dims[1] == 1){ bc[1] = 5;}
        else if (dims[1] > 1){
            if (cart_coords[1] == 0){            bc[1] = 1;}
            if (cart_coords[1] + 1 == dims[1]){  bc[1] = 3;}

        }

        if (dims[2] == 1){ bc[2] = 5;}
        else if (dims[2] > 1){
            if (cart_coords[2] == 0){            bc[2] = 1;}
            if (cart_coords[2] + 1 == dims[2]){  bc[2] = 3;}

        }
    }

    else if (bc_method == 4){
        if (dims[0] == 1){ bc[0] = 6;}
        else if (dims[0] > 1){
            if (cart_coords[0] == 0){            bc[0] = 2;}
            if (cart_coords[0] + 1 == dims[0]){  bc[0] = 4;}
        }

        if (dims[1] == 1){ bc[1] = 6;}
        else if (dims[1] > 1){
            if (cart_coords[1] == 0){            bc[1] = 2;}
            if (cart_coords[1] + 1 == dims[1]){  bc[1] = 4;}

        }

        if (dims[2] == 1){ bc[2] = 5;} // 5 here is right !
        else if (dims[2] > 1){
            if (cart_coords[2] == 0){            bc[2] = 2;}
            if (cart_coords[2] + 1 == dims[2]){  bc[2] = 4;}

        }
    }


}   


void set_halos_initial( double*** var, 
                        int* new_proc_size, 
                        unsigned int* bc)
{
    /*
    initially set halos to zero
    */
    int i, j, k;

    for (i = 0; i < new_proc_size[2]; i++){
        for (j = 0; j < new_proc_size[1]; j++){
            if (bc[0] == 2 || bc[0] == 6){
                var[i][j][0] = 0.0;
                var[i][j][0] = 0.0;
            }
            if (bc[0] == 4 || bc[0] == 6){
                var[i][j][new_proc_size[0]-1] = 0.0;
                var[i][j][new_proc_size[0]-2] = 0.0;
            }
        }
    }

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

    for (i = 0; i < new_proc_size[1]; i++){
        for (j = 0; j < new_proc_size[0]; j++){
            if (bc[3] == 2 || bc[3] == 6){
                var[0][i][j] = 0.0;
                var[1][i][j] = 0.0;
            }
            if (bc[3] == 4 || bc[3] == 6){
                var[new_proc_size[2]-1][i][j] = 0.0;
                var[new_proc_size[2]-2][i][j] = 0.0;
            }
        }
    }
}


