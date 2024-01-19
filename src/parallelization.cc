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
 * Purpose:    void communicate_geom_halos:
 *             used once, communicate halos geometry
 *
 *             void communicate_halos:
 *             pressure and velocity
 *             all 3 directions, if no communication is 
 *             needed due to BCs or small number of ranks
 *             MPI_PROC_NULL means that there is no neighbour 
 *             (defined as -2 in MPI4)
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0

              krach 2021/01/14 added different cases for bc
                                added pindicator since pressure-bc have 
                                to be treated different than velocity bc 
                                in case of no slip conditions
                                if pindicator == 0 --> var is pressure
                                if pindicator == 1 --> var is velocity
 ***********************************************************************/


// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "constants.h"

// **********************************************************************

void communicate_geom_halos(bool*** var, 
                            int* new_proc_size, 
                            int* neighbors, 
                            MPI_Datatype etype, 
                            MPI_Comm comm_cart)
{
    int         i, j, k;
    MPI_Status  srstatus[12];           // MPI_Status object for all communication dirs
    MPI_Datatype zsend, ysend, xsend;   // derived datatypes for communication

    int zcount = new_proc_size[0];
    int ycount = new_proc_size[0];
    int xcount = 1;

    MPI_Type_contiguous(zcount, etype, &zsend);
    MPI_Type_commit(&zsend);
    MPI_Type_contiguous(ycount, etype, &ysend);
    MPI_Type_commit(&ysend);
    MPI_Type_contiguous(xcount, etype, &xsend);
    MPI_Type_commit(&xsend);

    // comunicate in x direction <--> size[0], neighbor[0]: -1 in xdir, neighbor[1]: +1 in xdir
    // communication bool by bool is improveable
    for (i = 0; i < new_proc_size[2]; i++){
        for (j = 0; j < new_proc_size[1]; j++){
            MPI_Sendrecv(&var[i][j][2],                  1, xsend, neighbors[0], 0, 
                         &var[i][j][new_proc_size[0]-2], 1, xsend, neighbors[1], 0, comm_cart, &srstatus[0]);
            MPI_Sendrecv(&var[i][j][3],                  1, xsend, neighbors[0], 1, 
                         &var[i][j][new_proc_size[0]-1], 1, xsend, neighbors[1], 1, comm_cart, &srstatus[1]);
            MPI_Sendrecv(&var[i][j][new_proc_size[0]-4], 1, xsend, neighbors[1], 2, 
                         &var[i][j][0],                  1, xsend, neighbors[0], 2, comm_cart, &srstatus[2]);
            MPI_Sendrecv(&var[i][j][new_proc_size[0]-3], 1, xsend, neighbors[1], 3, 
                         &var[i][j][1],                  1, xsend, neighbors[0], 3, comm_cart, &srstatus[3]);
        }
    }

    // communicate in y direction <--> size[1], neighbor[2]: -1 in ydir, neighbor[3]: +1 in ydir
    for (i = 0; i < new_proc_size[2]; i++){
        // synchronous blocking sendrecvs
        MPI_Sendrecv(var[i][2],                  1, ysend, neighbors[2], 4,
                     var[i][new_proc_size[1]-2], 1, ysend, neighbors[3], 4, comm_cart, &srstatus[4]);
        MPI_Sendrecv(var[i][3],                  1, ysend, neighbors[2], 5,
                     var[i][new_proc_size[1]-1], 1, ysend, neighbors[3], 5, comm_cart, &srstatus[5]);
        MPI_Sendrecv(var[i][new_proc_size[1]-4], 1, ysend, neighbors[3], 6,
                     var[i][0],                  1, ysend, neighbors[2], 6, comm_cart, &srstatus[6]);
        MPI_Sendrecv(var[i][new_proc_size[1]-3], 1, ysend, neighbors[3], 7,
                     var[i][1],                  1, ysend, neighbors[2], 7, comm_cart, &srstatus[7]);
    }

    // communicate in z direction <--> size[2], neighbor[4]: -1 in zdir, neighbor[5]: +1 in zdir
    for (i = 0; i < new_proc_size[1]; i++){
        // synchronous blocking sendrecvs
        MPI_Sendrecv(var[2][i],                  1, zsend, neighbors[4], 8,
                     var[new_proc_size[2]-2][i], 1, zsend, neighbors[5], 8,  comm_cart, &srstatus[8]);
        MPI_Sendrecv(var[3][i],                  1, zsend, neighbors[4], 9, 
                     var[new_proc_size[2]-1][i], 1, zsend, neighbors[5], 9,  comm_cart, &srstatus[9]);
        MPI_Sendrecv(var[new_proc_size[2]-4][i], 1, zsend, neighbors[5], 10,
                     var[0][i],                  1, zsend, neighbors[4], 10, comm_cart, &srstatus[10]);
        MPI_Sendrecv(var[new_proc_size[2]-3][i], 1, zsend, neighbors[5], 11,
                     var[1][i],                  1, zsend, neighbors[4], 11, comm_cart, &srstatus[11]);
    }

}


void communicate_halos( double*** var, 
                        int* new_proc_size, 
                        int* neighbors, 
                        MPI_Datatype etype, 
                        MPI_Comm comm_cart, 
                        unsigned int* bc, 
                        int pindicator)
{
    int         i, j, k;
    MPI_Status  srstatus[12];           // MPI_Status object for all communication dirs
    MPI_Datatype zsend, ysend, xsend;   // derived datatypes for communication

    int zcount = new_proc_size[0];
    int ycount = new_proc_size[0];
    int xcount = 1;
    // int zmax_procs = dims[2] - 1; // max number of process in z dir in cart communicator

    // Definition of Datatypes to send/recv
    MPI_Type_contiguous(zcount, etype, &zsend);
    MPI_Type_commit(&zsend);
    MPI_Type_contiguous(ycount, etype, &ysend);
    MPI_Type_commit(&ysend);
    MPI_Type_contiguous(xcount, etype, &xsend);
    MPI_Type_commit(&xsend);

    // comunicate in x direction <--> size[0], neighbor[0]: -1 in xdir, neighbor[1]: +1 in xdir
    // communication bool by bool is improveable
    for (i = 0; i < new_proc_size[2]; i++){
        for (j = 0; j < new_proc_size[1]; j++){
            MPI_Sendrecv(&var[i][j][2],                  1, xsend, neighbors[0], 0, 
                         &var[i][j][new_proc_size[0]-2], 1, xsend, neighbors[1], 0, comm_cart, &srstatus[0]);
            MPI_Sendrecv(&var[i][j][3],                  1, xsend, neighbors[0], 1, 
                         &var[i][j][new_proc_size[0]-1], 1, xsend, neighbors[1], 1, comm_cart, &srstatus[1]);
            MPI_Sendrecv(&var[i][j][new_proc_size[0]-4], 1, xsend, neighbors[1], 2, 
                         &var[i][j][0],                  1, xsend, neighbors[0], 2, comm_cart, &srstatus[2]);
            MPI_Sendrecv(&var[i][j][new_proc_size[0]-3], 1, xsend, neighbors[1], 3, 
                         &var[i][j][1],                  1, xsend, neighbors[0], 3, comm_cart, &srstatus[3]);
        }
    }

    // communicate in y direction <--> size[1], neighbor[2]: -1 in ydir, neighbor[3]: +1 in ydir
    for (i = 0; i < new_proc_size[2]; i++){
        // synchronous blocking sendrecvs
        MPI_Sendrecv(var[i][2],                  1, ysend, neighbors[2], 4,
                     var[i][new_proc_size[1]-2], 1, ysend, neighbors[3], 4, comm_cart, &srstatus[4]);
        MPI_Sendrecv(var[i][3],                  1, ysend, neighbors[2], 5,
                     var[i][new_proc_size[1]-1], 1, ysend, neighbors[3], 5, comm_cart, &srstatus[5]);
        MPI_Sendrecv(var[i][new_proc_size[1]-4], 1, ysend, neighbors[3], 6,
                     var[i][0],                  1, ysend, neighbors[2], 6, comm_cart, &srstatus[6]);
        MPI_Sendrecv(var[i][new_proc_size[1]-3], 1, ysend, neighbors[3], 7,
                     var[i][1],                  1, ysend, neighbors[2], 7, comm_cart, &srstatus[7]);
    }

    // communicate in z direction <--> size[2], neighbor[4]: -1 in zdir, neighbor[5]: +1 in zdir
    for (i = 0; i < new_proc_size[1]; i++){

        // synchronous blocking sendrecvs
        MPI_Sendrecv(var[2][i],                  1, zsend, neighbors[4], 8,
                     var[new_proc_size[2]-2][i], 1, zsend, neighbors[5], 8,  comm_cart, &srstatus[8]);
        MPI_Sendrecv(var[3][i],                  1, zsend, neighbors[4], 9, 
                     var[new_proc_size[2]-1][i], 1, zsend, neighbors[5], 9,  comm_cart, &srstatus[9]);
        MPI_Sendrecv(var[new_proc_size[2]-4][i], 1, zsend, neighbors[5], 10,
                     var[0][i],                  1, zsend, neighbors[4], 10, comm_cart, &srstatus[10]);
        MPI_Sendrecv(var[new_proc_size[2]-3][i], 1, zsend, neighbors[5], 11,
                     var[1][i],                  1, zsend, neighbors[4], 11, comm_cart, &srstatus[11]);
        
    }

    // Free the datatypes 
    MPI_Type_free(&zsend);
    MPI_Type_free(&ysend);
    MPI_Type_free(&xsend);

}

