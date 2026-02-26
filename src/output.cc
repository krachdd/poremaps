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
 * Date:       2021/12
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    void write_output_raw:
 *             overloaded, analog to read file in parallel
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
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
#include <algorithm>
#include "mpi.h"
#include "constants.h"
// **********************************************************************

void write_output_raw(  double*** var, 
                        char* file_name, 
                        int* size, 
                        int* proc_size, 
                        int* new_proc_size, 
                        int* starts, 
                        MPI_Comm comm_cart, 
                        int* dims, 
                        int rank, 
                        int n_proc, 
                        int* cart_coords, 
                        MPI_Datatype etype, 
                        bool rescale, 
                        double voxelsize)
{

    int i, j, k;
    MPI_File      fh;  // File handle for var File
    int           order = MPI_ORDER_FORTRAN;
    MPI_Datatype  filetype;
    MPI_Offset    write_offset = 0;
    MPI_Status    write_status;

    MPI_Type_create_subarray(3, size, proc_size, starts, order, etype, &filetype);
    MPI_Type_commit(&filetype);
    // Open file in each rank of Cartesian Communicator (must be an intracomm)
    MPI_File_open(comm_cart, file_name, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    double *temp_var = new double[proc_size[2]*proc_size[1]*proc_size[0]];
    long int counter = 0;
    for (i = 0; i < proc_size[2]; i++){
        for (j = 0; j < proc_size[1]; j++){
            for (k = 0; k < proc_size[0]; k++){
                if (rescale == true){temp_var[counter] = var[i+2][j+2][k+2] * voxelsize;}
                else {temp_var[counter] = var[i+2][j+2][k+2];}
                counter++;
            }
        }
    }
    
    MPI_File_set_view(fh, write_offset, etype, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_at_all(  fh,
                            write_offset,
                            temp_var,
                            proc_size[2]*proc_size[1]*proc_size[0],
                            etype,
                            &write_status);

    delete [] temp_var;
    MPI_File_close(&fh);
    MPI_Type_free(&filetype);
}

void write_output_raw(  int*** var,
                        char* file_name, 
                        int* size,
                        int* proc_size,
                        int* new_proc_size,
                        int* starts,
                        MPI_Comm comm_cart,
                        int* dims,
                        int rank,
                        int n_proc,
                        int* cart_coords,
                        MPI_Datatype etype)
{

    int i, j, k;
    MPI_File      fh;  // File handle for var File
    int           order = MPI_ORDER_FORTRAN;
    MPI_Datatype  filetype;
    MPI_Offset    write_offset = 0;
    MPI_Status    write_status;

    MPI_Type_create_subarray(3, size, proc_size, starts, order, etype, &filetype);
    MPI_Type_commit(&filetype);
    // Open file in each rank of Cartesian Communicator (must be an intracomm)
    MPI_File_open(comm_cart, file_name, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    int *temp_var = new int[proc_size[2]*proc_size[1]*proc_size[0]];
    long int counter = 0;
    for (i = 0; i < proc_size[2]; i++){
        for (j = 0; j < proc_size[1]; j++){
            for (k = 0; k < proc_size[0]; k++){
                temp_var[counter] = var[i+2][j+2][k+2];
                counter++;
            }
        }
    }
    MPI_File_set_view(fh, write_offset, etype, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_at_all(  fh,
                            write_offset,
                            temp_var,
                            proc_size[2]*proc_size[1]*proc_size[0],
                            etype,
                            &write_status);

    delete [] temp_var;
    MPI_File_close(&fh);
    MPI_Type_free(&filetype);
}


void write_domain_decomposition(int*** domain_decomposition,
                                char* file_name,
                                int* size,
                                int* proc_size,
                                int* new_proc_size,
                                int* starts,
                                MPI_Comm comm_cart,
                                int* dims,
                                int rank,
                                int n_proc,
                                int* cart_coords,
                                MPI_Datatype etype)
{
    int i, j, k;
    for (i = 0; i < new_proc_size[2]; i++){
        domain_decomposition[i] = (int**)malloc((new_proc_size[1])*sizeof(int*));
        for (j = 0; j < new_proc_size[1]; j++){
            domain_decomposition[i][j] = (int*)malloc((new_proc_size[0])*sizeof(int));
            for (k = 0; k < new_proc_size[0]; k++){
                domain_decomposition[i][j][k] = rank;
            }
        }
    }
    write_output_raw(   domain_decomposition, 
                        file_name, 
                        size, 
                        proc_size, 
                        new_proc_size, 
                        starts, comm_cart, 
                        dims, 
                        rank, 
                        n_proc, 
                        cart_coords, 
                        etype);
}


void write_logfile( int iteration, 
                    double k13, 
                    double k23, 
                    double k33, 
                    double wk13, 
                    double wk23, 
                    double wk33,
                    double wmax_velz,
                    double wmean_velz, 
                    double conv, 
                    char* file_name, 
                    bool firstline, 
                    double tps)
{
    
    FILE* log_file;
    log_file = fopen(file_name, "a");
    


    if (firstline == true){
        fprintf(log_file, "# iteration, conv, TPS [1/s], wmax_velz [m/s], wmean_velz [m/s], k13 [m^2], k23 [m^2], k33 [m^2], wk13 [m^2], wk23 [m^2], wk33 [m^2]\n");
    }
    fprintf(log_file, "%i, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", 
                iteration,
                conv, 
                tps,
                wmax_velz,
                wmean_velz, 
                k13, 
                k23, 
                k33, 
                wk13, 
                wk23, 
                wk33 
                );
    fclose(log_file);
}
