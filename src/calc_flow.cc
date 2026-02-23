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
 * Purpose:    void calc_flow:
 *             main loop containing field updates, communication 
 *             and convergence check
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
#include "parallelization.h"
#include "solver.h"
#include "boundary_conditions.h"
#include "geometry.h"
#include "evaluation.h"
#include "output.h"
// **********************************************************************

void calc_flow(    bool*** proc_geom, 
                   int* size, 
                   int* new_proc_size, 
                   unsigned int* bc, 
                   double* dx, 
                   int max_iterations, 
                   double*** press, 
                   double*** vel_x, 
                   double*** vel_y, 
                   double*** vel_z,
                   int it_write, 
                   int it_eval, 
                   unsigned int solver,
                   double eps, 
                   double porosity,
                   int n_proc, 
                   int rank, 
                   int*** voxel_neighborhood, 
                   MPI_Comm comm_cart, 
                   int* dims, 
                   int* cart_coords, 
                   int* neighbors,
                   int* starts, 
                   int* ends,
                   int* dom_interest,
                   char* log_file)
{

    int iteration = 0;                          // current iteration counter
    double perm_conv = 1.0 + eps;               // convergence criterion (initialised above threshold)
    double start_it0, end_it0, elapsed_seconds; // per-iteration timing (rank 0 only)
    double phi_global;                          // domain porosity
    double permeability = 0.0;                  // mean z-velocity based permeability (voxel units)
    double wmax_velz = 0.0, wmean_velz = 0.0;   // whole-domain velocity diagnostics
    double k13 = 0.0, k23 = 0.0, k33 = 0.0;    // permeability tensor, domain of interest [m^2]
    double wk13 = 0.0, wk23 = 0.0, wk33 = 0.0; // permeability tensor, whole domain [m^2]
    bool   firstline = true;                    // controls whether to write the header line in the log file

    // if porosity is given as -1.0 --> compute it by get_porosity
    if (porosity == -1.0){
        phi_global = get_porosity(proc_geom, new_proc_size, comm_cart);
    }
    else { // Use the given value, useful if e.g. solid frame is added
        phi_global = porosity;
    }

    if (rank == 0){printf("START MAIN COMPUTATION\n");}

    while (iteration < max_iterations && perm_conv >= eps){
        // Main Loop for computation

        if (rank == 0){start_it0 = MPI_Wtime();}
        if (iteration == 0){
            // Set inital pressure gradient
            reapply_pressure_gradient(proc_geom, size, new_proc_size, press, comm_cart, cart_coords, dims);
        }
        central_diff5(proc_geom, new_proc_size, press, vel_x, vel_y, vel_z, voxel_neighborhood, 
                      solver);

        communicate_halos(press, new_proc_size, neighbors, MPI_DOUBLE, comm_cart, bc, 0);
        communicate_halos(vel_x, new_proc_size, neighbors, MPI_DOUBLE, comm_cart, bc, 1);
        communicate_halos(vel_y, new_proc_size, neighbors, MPI_DOUBLE, comm_cart, bc, 1);
        communicate_halos(vel_z, new_proc_size, neighbors, MPI_DOUBLE, comm_cart, bc, 1);
        // Set Dirichlet Boundarys for the pressure again
        reapply_pressure_gradient(proc_geom, size, new_proc_size, press, comm_cart, cart_coords, dims);

        // evaluate convergence
        if (iteration%it_eval == 0 && iteration > 1000){
            perm_conv = compute_convergence(proc_geom, 
                                            size, 
                                            new_proc_size, 
                                            vel_z, 
                                            &permeability, 
                                            comm_cart);

            // write output to logfile
            if ( iteration%it_write == 0 ){
                compute_permeability(   proc_geom,
                                            size,
                                            new_proc_size,
                                            vel_x, vel_y, vel_z,
                                            press, 
                                            &k13, &k23, &k33,
                                            &wk13, &wk23, &wk33,
                                            &wmax_velz, &wmean_velz,
                                            comm_cart, 
                                            starts,
                                            dom_interest,
                                            *dx,
                                            phi_global);

                if ( rank == 0 ){
                    end_it0 = MPI_Wtime(); elapsed_seconds = end_it0 - start_it0;
                    write_logfile(  iteration,
                                    k13, k23, k33,
                                    wk13, wk23, wk33,
                                    wmax_velz, wmean_velz,
                                    perm_conv,
                                    log_file,
                                    firstline, 
                                    (1./elapsed_seconds));
                    firstline = false;
                }
            }
        }

        if (rank == 0 && iteration%100 == 0){
            end_it0 = MPI_Wtime(); elapsed_seconds = end_it0 - start_it0;
            printf("Iteration: %i\t Convergence: %e\t", iteration, perm_conv);
            printf("Time: %f\t TPS: %f\n", elapsed_seconds, (1./elapsed_seconds));
        }
        // update iteration
        iteration++;
    

    }
} 

