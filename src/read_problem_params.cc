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
 * Date:       2021/02
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    void read_problem_params:
 *             Read problem parameters. File is always input_param.inp
 *             if not specified otherwise. 
 * 
 *             void read_parallel_params:
 *             Read parallel parameters, same function type, reads only
 *             parameters needed for domain decomposition
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 * 
 ***********************************************************************/

// HEADER ***************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits>
#include "mpi.h"
// **********************************************************************

void read_problem_params(   char* file_name, 
                            char* GEOM_FILE, 
                            int* size, 
                            double* dx, 
                            int* max_iterations, 
                            int* it_eval, 
                            int* it_write, 
                            char* log_file, 
                            unsigned int* solver, 
                            double* eps,
                            double* porosity,
                            int* dom_interest,
                            int* write_output,
                            int rank)
{   
    // Dummys for domain decomp. and bc_method to keep reading as simple as possible
    int  dummydims[3];
    int  dummybcs; 
    
    // Open file and set position indicator to beginning
    FILE* input_param_file;
    input_param_file = fopen(file_name, "r");
    rewind(input_param_file);

    // Read parameters
    fscanf(input_param_file,"dom_decomposition %d %d %d \n", &dummydims[0], &dummydims[1], &dummydims[2]);
    fscanf(input_param_file,"boundary_method %d \n", &dummybcs);
    fscanf(input_param_file,"geometry_file_name %s \n" ,GEOM_FILE);
    fscanf(input_param_file,"size_x_y_z %d %d %d \n" ,&size[0], &size[1], &size[2]); // global size 
    fscanf(input_param_file,"voxel_size  %lf \n", dx);
    fscanf(input_param_file,"max_iter %d \n", max_iterations);
    fscanf(input_param_file,"it_eval %d \n", it_eval);
    fscanf(input_param_file,"it_write %d \n", it_write);
    fscanf(input_param_file,"log_file_name %s \n", log_file);
    fscanf(input_param_file,"solving_algorithm %u \n", solver);
    fscanf(input_param_file,"eps %lf \n", eps);
    fscanf(input_param_file,"porosity %lf \n", porosity);
    fscanf(input_param_file,"dom_interest %d %d %d %d %d %d \n", &dom_interest[0], &dom_interest[1],
                                                                 &dom_interest[2], &dom_interest[3],
                                                                 &dom_interest[4], &dom_interest[5]);
    fscanf(input_param_file,"write_output %d %d %d %d \n", &write_output[0], &write_output[1],
                                                           &write_output[2], &write_output[3]);

    fclose(input_param_file); 

    // Sanity checks
    if (*eps <= 0){
        fprintf(stderr, "Setting EPS to a value smaller equal than zero.");
    }
    if (*it_write % *it_eval){
        fprintf(stderr, "Set writing and evaluating freq to proper values (it_write mod it_eval != 0).");
    }
    if (*porosity > 1.0 || (*porosity < -0.0 && *porosity != -1.0)){
        fprintf(stderr, "Porosity has to be 0.0 < porosity <= 1.0 or -1.0 to be evaluated.");
    }

    bool print_all_output = false;

    if (rank == 0){
        if (print_all_output == true ){
            printf("geometry_file_name %s\n", GEOM_FILE);
            printf("domain size %d, %d, %d\n", size[0], size[1], size[2]);
            printf("voxel size %lf\n", *dx);
            printf("maximum iterations %d\n", *max_iterations);
            printf("it_eval %d \n", *it_eval);
            printf("it_write %d \n", *it_write);
            printf("log_file_name %s \n", log_file);
            printf("solving_algorithm %u \n", *solver);
            printf("eps %lf \n", *eps);
            printf("porosity %lf \n", *porosity);
            printf("dom_interest %d %d %d %d %d %d \n", dom_interest[0], dom_interest[1],
                                                        dom_interest[2], dom_interest[3],
                                                        dom_interest[4], dom_interest[5]);
        }
    }
}


void read_parallel_params(  char* file_name, 
                            int* bc_method, 
                            int* dims,
                            int rank)
{   

    // Open file
    FILE* input_param_file;
    input_param_file = fopen(file_name, "r");

    fscanf(input_param_file,"dom_decomposition %d %d %d \n", &dims[0], &dims[1], &dims[2]);
    fscanf(input_param_file,"boundary_method %u \n", bc_method);

    fclose(input_param_file); 

    // Sanity checks
    if (!(dims[0] == 0 && dims[1] == 0 && dims[2] == 0) && !(dims[0] != 0 && dims[1] != 0 && dims[2] != 0)){
        fprintf(stderr, "Check domain decomposition.");
    }

    bool print_all_output = false;

    if (rank == 0){
        if (print_all_output == true ){
            printf("dom_decomposition %d %d %d \n", dims[0], dims[1], dims[2]);
            printf("boundary_method %u \n", *bc_method);
        }
    }
}