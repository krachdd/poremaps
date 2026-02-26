/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright 2024-2026 David Krach, Matthias Ruf

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
 * Date:       2021/02
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    Evaluate Geometry. Neighbourhood. Read input geometry file (.raw)
 * 
 *             void read_geometry: 
 *             read file in parallel at all ranks
 * 
 *             void get_dom_limits:
 *             compute the domain limits dependend on number
 *             of ranks and domain decomposition
 *
 * Contents:   Version: 1.0.0
 * 
 ***********************************************************************/

// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iterator>
#include <numeric>
#include <math.h>
#include <vector>

#include "mpi.h"
#include "constants.h"
#include "parallelization.h"
// **********************************************************************


void read_geometry( const char* GEOM_FILE, 
                    int* size, 
                    int* proc_size, 
                    int* new_proc_size, 
                    int* starts, 
                    bool*** proc_geom, 
                    MPI_Comm comm_cart, 
                    int* dims, 
                    int rank, 
                    int n_proc, 
                    int* cart_coords,
                    int* neighbors)
{
    int i, j, k;
    MPI_File      geom_fh;                      // File handle for Geometry File
    int           order = MPI_ORDER_FORTRAN;    // Order of file in storage
    MPI_Datatype  etype = MPI_C_BOOL;           // Datatype of file in storage
    MPI_Datatype  filetype;                     // Derived datatype
    MPI_Offset    read_offset;                  // Helps to get absolute byte position
    MPI_Status    read_status; 


    if (rank == 0){printf("Reading input geometry\n");}
    // Open file in each rank of Cartesian Communicator (must be an intracomm)
    MPI_Type_create_subarray(3, size, proc_size, starts, order, etype, &filetype);
    MPI_Type_commit(&filetype);
    // Do not open file with MPI_MODE_SEQUENTIAL
    MPI_File_open(comm_cart, GEOM_FILE, MPI_MODE_RDONLY, MPI_INFO_NULL, &geom_fh);
    read_offset = 0;
    MPI_File_set_view(geom_fh, read_offset, etype, filetype, "native", MPI_INFO_NULL);

    bool *temp_geom = new bool[proc_size[2]*proc_size[1]*proc_size[0]];

    MPI_File_read_at_all(   geom_fh, 
                            read_offset, 
                            temp_geom, 
                            proc_size[0]*proc_size[1]*proc_size[2], 
                            etype, 
                            &read_status);

    for (i = 0; i < new_proc_size[2]; i++){
        proc_geom[i] = (bool**)malloc((new_proc_size[1])*sizeof(bool*));
        for (j = 0; j < new_proc_size[1]; j++){
            proc_geom[i][j] = (bool*)malloc((new_proc_size[0])*sizeof(bool));
            for (k = 0; k < new_proc_size[0]; k++){
                proc_geom[i][j][k] = 0;
            }
        }
    }

    // long int to be absolutly sure
    long int counter = 0;
    // shift domain by halo to central location
    for (i = 0; i < proc_size[2]; i++){
        for (j = 0; j < proc_size[1]; j++){
            for (k = 0; k < proc_size[0]; k++){
                proc_geom[i+2][j+2][k+2] = temp_geom[counter];
                counter++;
            }
        }
    }

    // Communicate halos
    communicate_geom_halos(proc_geom, new_proc_size, neighbors, etype, comm_cart);
    MPI_File_close(&geom_fh);

    // free array on heap 
    delete [] temp_geom;
}


void get_dom_limits(    int* size, 
                        int* proc_size, 
                        int* new_proc_size, 
                        int* starts, 
                        int* ends, 
                        int* dims, 
                        int rank, 
                        int* cart_coords)
{
    int i; 
    int addslice[3]; // to deal with modulo 
    int cum_addslices; // cummulative number of to be added slices


    // determine problem size per process within Cartesian Topology
    for (i = 0; i < 3; i++){proc_size[i] = floor(size[i]/dims[i]);}
    for (i = 0; i < 3; i++){addslice[i] = size[i] - dims[i] * proc_size[i];}

    // compute domain limits
    for (i = 0; i < 3; i++){
        // compute addslices 
        if (cart_coords[i] == 0){
            cum_addslices = 0;
        }
        else if (cart_coords[i] > 0 && cart_coords[i] <= addslice[i]){
            cum_addslices = cart_coords[i];
        }
        else {
            cum_addslices = addslice[i];
        }
        starts[i] = cart_coords[i]*proc_size[i] + cum_addslices;
    }
    for (i = 0; i < 3; i++){ends[i] = starts[i] + proc_size[i];}
    // compute new problem sizes on procs
    for (i = 0; i < 3; i++){
        // if (rank < addslice[i]){
        if (cart_coords[i] < addslice[i]){
            proc_size[i] = proc_size[i] + 1;
        }
    }
    // Add halos
    for (i = 0; i < 3; i++){new_proc_size[i] = proc_size[i] + 4;}
}


void eval_geometry( bool*** proc_geom, 
                    int* new_proc_size, 
                    int*** voxel_neighborhood)
{
    // See Bentz et al 2007.
    int i, j, k, direction;
    int i_, j_, k_;
    int shift_factor = 0;


    // allocate memory for voxel_neighborhood with halos on each process
    for (i = 0; i < new_proc_size[2]; i++){
        voxel_neighborhood[i] = (int**)malloc((new_proc_size[1])*sizeof(int*));
        for (j = 0; j < new_proc_size[1]; j++){
            voxel_neighborhood[i][j] = (int*)malloc((new_proc_size[0])*sizeof(int));
            for (k = 0; k < new_proc_size[0]; k++){
                // set inital nhood to zero if voxel is fluid
                if (proc_geom[i][j][k] == 0){
                    voxel_neighborhood[i][j][k] = 0;
                }
                // set inital nhood to -1 if voxel is solid
                else {
                    voxel_neighborhood[i][j][k] = -1;
                }
            }
        }
    }

    for (direction = 0; direction < 3; direction++){
        // x direction
        if (direction == 0) {
            i_ = 0;
            j_ = 0;
            k_ = 1;
            shift_factor = 100;
        }
        // y direction
        else if (direction == 1) {
            i_ = 0;
            j_ = 1;
            k_ = 0;
            shift_factor = 10;
        }
        // z direction
        else if (direction == 2) {
            i_ = 1;
            j_ = 0;
            k_ = 0;
            shift_factor = 1;
        }

        int lim = 2;

        for (i = lim; i < new_proc_size[2]-lim; i++){
            for (j = lim ; j < new_proc_size[1]-lim; j++){
                for (k = lim; k < new_proc_size[0]-lim; k++){
                    // definitions used for voxel check
                    // fluid voxel :    =  0
                    // solid voxel :    =  1
                    // outside domain : =  2

                    int lvox = -1;      // left voxel
                    int llvox = -1;     // second left voxel
                    int rvox = -1;      // right voxel
                    int rrvox = -1;     // second right voxel

                    // check only fluid voxels
                    if (proc_geom[i][j][k] == 0){
                        if (proc_geom[i - i_][j - j_][k - k_] == 0) {
                            lvox = 0;
                        }
                        else {
                            lvox = 1;
                        }
                        if (proc_geom[i + i_][j + j_][k + k_] == 0) {
                            rvox = 0;
                        }
                        else {
                            rvox = 1;
                        }

                        if (lvox == 0){
                            if (proc_geom[i - 2*i_][j - 2*j_][k - 2*k_] == 0) {
                                llvox = 0;
                            }
                            else {
                                llvox = 1;
                            }
                        }

                        if (rvox == 0) {
                            if (proc_geom[i + 2*i_][j + 2*j_][k + 2*k_] == 0) {
                                rrvox = 0;
                            }
                            else {
                                rrvox = 1;
                            }
                        }

                        // case 1
                        if ((lvox == 1) && (rvox == 1)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 1 * shift_factor;
                        }
                        // case 2
                        else if ((lvox == 1) && (rvox == 0) && (rrvox == 1)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 2 * shift_factor;
                        }
                        // case 3
                        else if ((llvox == 1) && (lvox == 0) && (rvox == 1)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 3 * shift_factor;
                        }
                        // case 4
                        else if ((lvox == 1) && (rvox == 0) && (rrvox == 0)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 4 * shift_factor;
                        }
                        // case 5
                        else if ((llvox == 0) && (lvox == 0) && (rvox == 1)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 5 * shift_factor;
                        }
                        // case 6
                        else if ((lvox == 0) && (rvox == 0)) {
                            voxel_neighborhood[i][j][k] = voxel_neighborhood[i][j][k] + 6 * shift_factor;
                        }                        
                    }
                }
            }
        }
    }
}


/* 
One has to be carefull if solid frame is added 
*/
double get_proc_porosity(   bool*** proc_geom, 
                            int* new_proc_size)
{
    int i, j, k;
    int fluid_count = 0;
    double porosity;
    int lim = 2;
    double vox_count = (new_proc_size[0] - 2*lim) * (new_proc_size[1] - 2*lim) * (new_proc_size[2] - 2*lim);

    // loop through whole geom
    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                // distinguish between solid and fluid voxels
                if (proc_geom[i][j][k] == 0){
                    fluid_count++;
                }
            }
        }
    }
    porosity = fluid_count/vox_count;
    return porosity;
}

/* 
One has to be carefull if solid frame is added 
*/
double get_porosity(    bool*** proc_geom, 
                        int* new_proc_size, 
                        MPI_Comm comm_cart)
{
    int i, j, k;
    double fluid_count = 0.0;
    double g_porosity;
    int lim = 2;
    double vox_count = (new_proc_size[0] - 2*lim) * (new_proc_size[1] - 2*lim) * (new_proc_size[2] - 2*lim);

    // loop through whole geom
    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                // distinguish between solid and fluid voxels
                if (proc_geom[i][j][k] == 0){
                    fluid_count += 1.0;
                }
            }
        }
    }

    // Batch both SUM reductions into one call.
    double gbuf[2] = {fluid_count, vox_count};
    MPI_Allreduce(MPI_IN_PLACE, gbuf, 2, MPI_DOUBLE, MPI_SUM, comm_cart);
    fluid_count = gbuf[0]; vox_count = gbuf[1];

    g_porosity = fluid_count/vox_count;
    return g_porosity;
}