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
 * Purpose:    void compute_permeability:
 *             computes everything to be put in logfile
 *
 *             double compute_convergence:
 *             compute convergence criterion
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

void compute_permeability(  bool*** proc_geom, 
                            int* size, 
                            int* new_proc_size, 
                            double*** vel_x, 
                            double*** vel_y, 
                            double*** vel_z,
                            double*** press, 
                            double* k13,
                            double* k23,
                            double* k33,
                            double* wk13,
                            double* wk23,
                            double* wk33,
                            double* wmax_velz, 
                            double* wmean_velz, 
                            MPI_Comm comm_cart,
                            int* starts, 
                            int* dom_interest,
                            double voxelsize,
                            double porosity)
{
    int i, j, k;
    int gi, gj, gk;                                     // local coordinates
    int lim = 2;                                        // to consider halos
    double velx_sum = 0.0, vely_sum = 0.0, velz_sum = 0.0;  
    double wvelx_sum = 0.0, wvely_sum = 0.0, wvelz_sum = 0.0;
    double wvelz_max = 0.0;                              // store velocities
    double p0sum = 0.0, pendsum = 0.0;                   // stores pressuers
    double pgrad_dofi = 1.0;                             // pressure gradient of subdomain
    int p0flcount = 0, pendflcount = 0;                  // voxel counter
    int flcount = 0, wflcount = 0;                       // voxel counter
    bool eval_dofi = true;                               // If domain of interest is evaluated

    for ( i = lim; i < new_proc_size[2] - lim; i++ ){
        for ( j = lim; j < new_proc_size[1] - lim; j++ ){
            for ( k = lim; k < new_proc_size[0] - lim; k++ ){
                // distinguish between solid and fluid voxels
                if ( proc_geom[i][j][k] == 0 ){
                    wflcount += 1;
                    wvelx_sum += vel_x[i][j][k];
                    wvely_sum += vel_y[i][j][k];
                    wvelz_sum += vel_z[i][j][k];
                    if ( vel_z[i][j][k] > wvelz_max ){
                        wvelz_max = vel_z[i][j][k];
                    }
                }
            }
        }
    }
    // if true --> compute properties of subdomain 
    if (dom_interest[0] == 0 && dom_interest[1] == 0 && 
        dom_interest[2] == 0 && dom_interest[3] == 0 && 
        dom_interest[4] == 0 && dom_interest[5] == 0){
        flcount = 0;
        eval_dofi = false;
    }
    else {
        // Check wether the subdomain of the rank is within the  
        // the domain domain of interest
        for ( i = lim; i < new_proc_size[2] - lim; i++ ){
            for ( j = lim; j < new_proc_size[1] - lim; j++ ){
                for ( k = lim; k < new_proc_size[0] - lim; k++ ){
                    // get the global location of the voxel
                    gi = starts[2] + i - lim;
                    gj = starts[1] + j - lim;
                    gk = starts[0] + k - lim;
                    if ( gk < dom_interest[0] || gk > dom_interest[1] ||
                         gj < dom_interest[2] || gj > dom_interest[3] ||
                         gi < dom_interest[4] || gk > dom_interest[5]   ){
                        continue;
                    }
                    else {
                        if ( proc_geom[i][j][k] == 0 ){
                            flcount += 1;
                            velx_sum += vel_x[i][j][k];
                            vely_sum += vel_y[i][j][k];
                            velz_sum += vel_z[i][j][k];
                        }
                    }
                    // To compute pressure gradient for the domain of interest
                    if ( gi == dom_interest[4] ){
                        if ( proc_geom[i][j][k] == 0 ){
                            p0flcount += 1;
                            p0sum += press[i][j][k];
                        }
                    }
                    if ( gi == dom_interest[5] ){
                        if ( proc_geom[i][j][k] == 0 ){
                            pendflcount += 1;
                            pendsum += press[i][j][k];
                        }
                    }
                }
            }
        }
    }

    // Make sure that there is no random value in storage and set all
    // relevant vels to zero if there is no fluid voxel present on the 
    // subdomain of the rank
    if ( wflcount == 0 ){
        wvelx_sum = 0.0;
        wvely_sum = 0.0;
        wvelz_sum = 0.0;
        wvelz_max = 0.0;
    }
    if ( flcount == 0 ){
        velx_sum = 0.0;
        vely_sum = 0.0;
        velz_sum = 0.0;
    }

    if ( p0flcount == 0 ){
        p0sum  = 0.0;
    }
    if ( pendflcount == 0 ){
        pendsum = 0.0;
    }

    // Communicate domain of interest properties if necessary 
    MPI_Barrier(comm_cart);

    // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
    if ( eval_dofi == true ){
        MPI_Allreduce(MPI_IN_PLACE, &velx_sum,    1, MPI_DOUBLE, MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &vely_sum,    1, MPI_DOUBLE, MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &velz_sum,    1, MPI_DOUBLE, MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &p0sum,       1, MPI_DOUBLE, MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &pendsum,     1, MPI_DOUBLE, MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &flcount,     1, MPI_INT,    MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &p0flcount,   1, MPI_INT,    MPI_SUM, comm_cart);
        MPI_Allreduce(MPI_IN_PLACE, &pendflcount, 1, MPI_INT,    MPI_SUM, comm_cart);
    }

    MPI_Barrier(comm_cart);

    MPI_Allreduce(MPI_IN_PLACE, &wvelx_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE, &wvely_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE, &wvelz_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE, &wflcount,  1, MPI_INT,    MPI_SUM, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE, &wvelz_max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);

    MPI_Barrier(comm_cart);

    // Compute permeabilities of domain of interest in real dimensions
    if ( eval_dofi == true ){
        pgrad_dofi = ((pendsum/pendflcount) - (p0sum/p0flcount))/(dom_interest[5] - dom_interest[4]);
        *k13 = -(velx_sum/flcount) * (1.0/Re) * porosity * (1.0/pgrad_dofi) * pow(voxelsize,2);
        *k23 = -(vely_sum/flcount) * (1.0/Re) * porosity * (1.0/pgrad_dofi) * pow(voxelsize,2);
        *k33 = (velz_sum/flcount) * (1.0/Re) * porosity * (1.0/pgrad_dofi) * pow(voxelsize,2);
    }
    else {
        *k13 = 0.0;
        *k23 = 0.0;
        *k33 = 0.0;
    }

    // Compute permeabilities and velocities of whole domain in physical dimensions
    *wk13 = -(wvelx_sum/wflcount) * (1.0/Re) * porosity * pow(voxelsize,2);
    *wk23 = -(wvely_sum/wflcount) * (1.0/Re) * porosity * pow(voxelsize,2);
    *wk33 =  (wvelz_sum/wflcount) * (1.0/Re) * porosity * pow(voxelsize,2);
    
    *wmean_velz = wvelz_sum/wflcount * voxelsize;
    *wmax_velz = wvelz_max * voxelsize;

    MPI_Barrier(comm_cart);

}




double compute_convergence( bool*** proc_geom, 
                            int* size, 
                            int* new_proc_size, 
                            double*** vel_z, 
                            double* permeability,
                            MPI_Comm comm_cart)
{
    int i, j, k;
    int lim = 2;                        // halo size to be excluded
    double mean_velz = 0.0;             // store mean vel in main dir
    double current_permeability = 0.0;  // permeability at current timestep
    double velz_sum = 0.0;              // sum of vel z of all fluid voxels
    double conv = 0.0;                  // convergence criterion
    int flcount = 0;                    // number of fluid voxels on rank

    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                // distinguish between solid and fluid voxels
                if (proc_geom[i][j][k] == 0){
                    flcount += 1;
                    velz_sum += vel_z[i][j][k];
                }
            }
        }
    }



    // Make sure that there is no random value in storage and set all
    // relevant vels to zero if there is no fluid voxel present on the 
    // subdomain of the rank
    if (flcount == 0){
        velz_sum = 0.0;
    }

    MPI_Barrier(comm_cart);

    // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
    MPI_Allreduce(MPI_IN_PLACE, &velz_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE, &flcount,  1, MPI_INT,    MPI_SUM, comm_cart);

    MPI_Barrier(comm_cart);

    mean_velz = velz_sum/flcount;
    current_permeability = mean_velz/Re; // in vx^2;
    // Compute convergence criterion
    conv = fabs(current_permeability - *permeability)/current_permeability;
    *permeability = current_permeability;
    MPI_Barrier(comm_cart);

    return conv;

}