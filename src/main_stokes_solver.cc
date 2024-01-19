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
 * Purpose:    Main Stokes-Solver File
 *             The programm is designed to analyze fluid flow patterns, 
 *             pressure fields and permeability in a 3D porous domain under
 *             a Stokes regime. Therefore a Finite Difference discretization
 *             is implemented.
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0 (for publication)
 ***********************************************************************/


// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "mpi.h"
#include "read_problem_params.h"
#include "constants.h"
#include "geometry.h"
#include "output.h"
#include "calc_flow.h"
#include "parallelization.h"
#include "initialize.h"
//**********************************************************************

int main (int argc, char** argv){
    /***********************************************************************
     * Declare all variables
    ************************************************************************/
    // parallelization
    int           rank, n_proc;       // Number of rank, number of processes
    int           rank_glob;          // Global number of rank, will be rearanged if reorder = true
    MPI_Comm      comm_cart;          // Communicator handle within Cartesian topology

    int           dims[ndims];        // specifying the number of processes in each dimension
    int           periods[ndims];     // Logical array of size ndims specifying whether the
                                      // grid is periodic (true) or not (false) in each dimension
    const int     reorder = true;     // Logical reordering of ranks in subcomm
                                      // If reorder = false then the rank of each process in the 
                                      // new group is identical to its rank in the old group. 
                                      // Otherwise, the function may reorder --> choose a good 
                                      // embedding of the virtual topology onto the physical machine
    int           neighbors[2*ndims]; // Array of ints to store nearest neighbor ranks (left, right, 
                                      // top, bottom, front, back)
    int           cart_coords[ndims]; // Cartesian coordinates of specified process

    // Parametrization 
    int           i;                  // 
    double        eps;                // Convergence Criterion 
    int           max_iterations;     // Breakup Criterion - maximal number of iterations 
    int           it_eval, it_write;  // frequency for evaluating progress and writing log file 
                                      // it_write mod it_eval =! 0 has to hold
    unsigned int  solver;             // Solver to be used (more information in header file)

    // File names
    std::string   input_fn;
    char          input_filename[file_name_char];
    char          GEOM_FILE[file_name_char];
    char          press_filename[file_name_char] = "press_";
    char          velx_filename[file_name_char] = "velx_";
    char          vely_filename[file_name_char] = "vely_";
    char          velz_filename[file_name_char] = "velz_";
    char          domd_filename[file_name_char] = "domain_decomp_";
    char          voxel_neighborhood_filename[file_name_char] = "voxel_neighborhood_";
    char          log_file[file_name_char];
    int           write_output[4];      // flag for writing output for vels, press, nhood, dom_decomp
                                        // if 1 write specified field to raw file

    // Problem variables 
    double        dx;                       // voxel size [m]
    double        porosity;                 // Porosity of the sample
    int           size[3];                  // Problem size (+ mirror domain + halo size)
    bool          ***proc_geometry;         // Geometry pointer
    int           ***voxel_neighborhood;    // Fluid voxel neighborhood pointer
    int           ***domain_decomposition;  // Domain Decomposition pointer
    unsigned int  bc[3];                    // Save BC per Direction (Documentation see initialize.cc)
    int           bc_method;                // get BC from Input file (Documentation see initialize.cc)
    int           proc_size[3];             // Subdomain Size per rank (halo etc. not included)
    int           new_proc_size[3];         // Real subdomain size per rank
    int           starts[3], ends[3];       // Starting/Ending coordinates of the subarray in each dimension 
                                            // (see MPI_Type_create_subarray documentation)
    int           dom_interest[2*ndims];    // xmin, xmax, ymin, ymax, zmin, zmax of subdomain for which
                                            // the permeability should be computed

    // Field variables
    double        ***press;           // Pressure
    double        ***vel_x;           // Velocity in x direction
    double        ***vel_y;           // Velocity in y direction
    double        ***vel_z;           // Velocity in z direction
    


    /***********************************************************************
     * End of Declare all variables
    ************************************************************************/
    
    // initialize MPI and check for number of processors
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_glob);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    if (rank_glob == 0){
        std::time_t starttime = std::time(nullptr);
        std::cout << "\nStarting Stokes Solver at " << std::asctime(std::localtime(&starttime));

        std::cout << "\nFinite Difference Stokes Solver  Copyright (C) 2023  David Krach, Matthias Ruf" << std::endl;
        std::cout << "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.md." << std::endl;
        std::cout << "This is free software, and you are welcome to redistribute it" << std::endl;
        std::cout << "under certain conditions; see LICENSE.md for details.\n" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Define input file
    // If no command line argument is given the solver assumes the input file is called "input_param.inp"
    if (argc == 1){input_fn = "input_param.inp";}
    else if (argc = 2){input_fn = (argv[1]);}

    if (rank_glob == 0){
        if (argc == 1){printf("No input file specified. Using input_param.inp\n");}
        else if (argc == 2){printf("Input file: %s\n", input_fn.c_str());}

        // Test if input file exists
        if (FILE *file = fopen(input_fn.c_str(), "r")) {fclose(file);}
        else {fprintf(stderr, "Input file doesnt exist.");}
    }

    memcpy(input_filename, input_fn.c_str(), input_fn.size()+1);


    // the read problems is necessary here since we need the bc_method to determine the periods
    read_parallel_params(input_filename, &bc_method, dims, rank);
    if (bc_method == 0){                  periods[0] = 1; periods[1] = 1; periods[2] = 1;}
    if (bc_method == 1 || bc_method == 2){periods[0] = 0; periods[1] = 0; periods[2] = 1;}
    if (bc_method == 3 || bc_method == 4){periods[0] = 0; periods[1] = 0; periods[2] = 0;}

    if ( dims[0] == 0 && dims[1] == 0 && dims[2] == 0){
        // if input is zero MPI determines domain decomposition
        for(i = 0; i < ndims; i++){dims[i] = 0;} // fill dims array
        MPI_Dims_create(n_proc, ndims, dims);
        if (rank_glob == 0){
            printf("MPI domain decomposition nx = %i, ny = %i, nz = %i\n", dims[0], dims[1], dims[2]);
        }
    }
    else {
        if (rank_glob == 0){
            printf("Given domain decomposition nx = %i, ny = %i, nz = %i\n", dims[0], dims[1], dims[2]);
        }
    }
    
    // Create Cartesian Topology within new MPI Communicator
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);
    // Get reordered ranks
    MPI_Comm_rank(comm_cart, &rank);
    // process coords in Cartesian topology given rank
    MPI_Cart_coords(comm_cart, rank, ndims, cart_coords);

    printf("Rank %i reports --> Coordinates: cx = %i, cy = %i, cz = %i\n", 
           rank, cart_coords[0], cart_coords[1], cart_coords[2]); 
    // Get shifted source and destination ranks 
    for(i = 0; i < ndims; i++){
        MPI_Cart_shift(comm_cart, i, 1, &neighbors[i*2], &neighbors[i*2+1]);
    }
    // Read problem parameters on each rank
    read_problem_params(    input_filename, 
                            GEOM_FILE, 
                            size, 
                            &dx, 
                            &max_iterations, 
                            &it_eval, 
                            &it_write, 
                            log_file, 
                            &solver, 
                            &eps,
                            &porosity,
                            dom_interest,
                            write_output,
                            rank);

    // Compute domain limits
    get_dom_limits(size, proc_size, new_proc_size, starts, ends, dims, rank, cart_coords);
    // Determine, depending on BC, which rank is and which is not communicating with whom
    determine_comm_pattern(dims, rank, cart_coords, bc, bc_method, neighbors);
    // MPI_PROC_NULL which means that there is no neighbour is defined as -2 in MPI4


    printf("Rank %i: Coordinates: cx = %i, cy = %i, cz = %i\n", 
                rank, cart_coords[0], cart_coords[1], cart_coords[2]);
    printf("Rank %i: Domain size: domx = %i, domy = %i, domz = %i\n", 
                rank, proc_size[0], proc_size[1], proc_size[2]);
    printf("Rank %i: Full domain size incl. halos: domx = %i, domy = %i, domz = %i\n", 
                rank, new_proc_size[0], new_proc_size[1], new_proc_size[2]);
    printf("Rank %i: Domain limits: x_lims (%i, %i), y_lims (%i, %i), z_lims (%i, %i)\n", 
                rank, starts[0], ends[0], starts[1], ends[1], starts[2], ends[2]);
    printf("Rank %i: Neighbors: left = %i, right = %i, top = %i, bottom = %i, front = %i, back = %i\n", 
                rank, neighbors[0], neighbors[1], neighbors[2], neighbors[3], neighbors[4], neighbors[5]);
    printf("Rank %i: Boundary Conditions: x : %i, y : %i, z : %i\n", 
                rank, bc[0], bc[1], bc[2]);
    MPI_Barrier(comm_cart);

    // allocate arrays
    proc_geometry      = (bool***)malloc((new_proc_size[2])*sizeof(bool**));
    voxel_neighborhood = (int***)malloc((new_proc_size[2])*sizeof(int**));
    press              = (double***)malloc((new_proc_size[2])*sizeof(double**));
    vel_x              = (double***)malloc((new_proc_size[2])*sizeof(double**));
    vel_y              = (double***)malloc((new_proc_size[2])*sizeof(double**));
    vel_z              = (double***)malloc((new_proc_size[2])*sizeof(double**));
    if ( write_output[3] == 1 ){
        domain_decomposition = (int***)malloc((new_proc_size[2])*sizeof(int**));
    }

    // Open input geometry 
    // Open file in each rank of Cartesian Communicator (must be an intracomm)
    read_geometry(  GEOM_FILE, 
                    size, 
                    proc_size, 
                    new_proc_size, 
                    starts, 
                    proc_geometry, 
                    comm_cart, 
                    dims, 
                    rank, 
                    n_proc, 
                    cart_coords, 
                    neighbors);

    MPI_Barrier(comm_cart);
    
    // Evaluate fluid voxel neighborhood - done only once
    if ( rank == 0 ){printf("Evaluate fluid voxel neighborhood\n");}
    eval_geometry(proc_geometry, new_proc_size, voxel_neighborhood);

    // Allocate arrays and set initial conditions
    allocate_initial_field( proc_geometry, 
                            size, 
                            new_proc_size, 
                            press, 
                            vel_x, vel_y, vel_z, 
                            comm_cart, 
                            neighbors, 
                            starts, 
                            bc);
    set_halos_initial(vel_x, new_proc_size, bc); 
    set_halos_initial(vel_y, new_proc_size, bc); 
    set_halos_initial(vel_z, new_proc_size, bc); 
    
    // main flow computation
    calc_flow(  proc_geometry, 
                size, 
                new_proc_size, 
                bc, 
                &dx, 
                max_iterations, 
                press, 
                vel_x, 
                vel_y, 
                vel_z, 
                it_write, 
                it_eval, 
                solver, 
                eps,
                porosity, 
                n_proc, 
                rank, 
                voxel_neighborhood, 
                comm_cart, 
                dims, 
                cart_coords, 
                neighbors,
                starts,
                ends,
                dom_interest, 
                log_file);

    // concatenate output names 
    strcat(voxel_neighborhood_filename, GEOM_FILE); 
    strcat(press_filename, GEOM_FILE); 
    strcat(velx_filename, GEOM_FILE); 
    strcat(vely_filename, GEOM_FILE); 
    strcat(velz_filename, GEOM_FILE);
    strcat(domd_filename, GEOM_FILE);
    

    if (rank == 0){printf("Writing output to file\n");}
    
    if ( write_output[2] == 1 ){
        //write results - voxel neighborhood
        if (rank == 0){printf("Writing neighborhood to file\n");}
        write_output_raw(   voxel_neighborhood, 
                            voxel_neighborhood_filename, 
                            size, 
                            proc_size, 
                            new_proc_size, 
                            starts, 
                            comm_cart, 
                            dims, 
                            rank, 
                            n_proc, 
                            cart_coords, 
                            MPI_INT);
    }
    
    if ( write_output[1] == 1 ){
        if (rank == 0){printf("Writing press to file\n");}
        write_output_raw(   press,
                            press_filename,
                            size, 
                            proc_size, 
                            new_proc_size, 
                            starts, 
                            comm_cart, 
                            dims, 
                            rank, 
                            n_proc, 
                            cart_coords, 
                            MPI_DOUBLE, 
                            false, 
                            dx);
    }

    if ( write_output[0] == 1 ){
        if (rank == 0){printf("Writing vels to file\n");}
        write_output_raw(   vel_x,
                            velx_filename, 
                            size, 
                            proc_size, 
                            new_proc_size, 
                            starts, 
                            comm_cart, 
                            dims, 
                            rank, 
                            n_proc, 
                            cart_coords, 
                            MPI_DOUBLE, 
                            true, 
                            dx);
        write_output_raw(   vel_y,
                            vely_filename, 
                            size, 
                            proc_size, 
                            new_proc_size, 
                            starts, 
                            comm_cart, 
                            dims, 
                            rank, 
                            n_proc, 
                            cart_coords, 
                            MPI_DOUBLE, 
                            true, 
                            dx);
        write_output_raw(   vel_z,
                            velz_filename, 
                            size, 
                            proc_size, 
                            new_proc_size, 
                            starts, 
                            comm_cart, 
                            dims, 
                            rank, 
                            n_proc, 
                            cart_coords, 
                            MPI_DOUBLE, 
                            true, 
                            dx);
    }
    if ( write_output[3] == 1 ){
        if ( rank == 0 ){printf("Writing domain decomposition to file\n");}
        // etype = MPI_INT;
        write_domain_decomposition( domain_decomposition, 
                                    domd_filename, 
                                    size, 
                                    proc_size, 
                                    new_proc_size, 
                                    starts, 
                                    comm_cart, 
                                    dims, 
                                    rank, 
                                    n_proc, 
                                    cart_coords, 
                                    MPI_INT);
    }

    // End MPI Processes

    MPI_Finalize();
} 
