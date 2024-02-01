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
 * Purpose:    void central_diff5:
 *             main part of pressure and vel update
 *              
 *             void get_2nd_derivation:
 *             2nd order derivatives based on neighborhood
 * 
 *             void split_number:
 *             helper function for neighborhood 
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0

 Solver:    1: Jacobi for the pressure
            2: Gauss-Seidel for pressure and fluid velocity
            3: SOR
 ***********************************************************************/


// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "constants.h"
// **********************************************************************


void get_2nd_derivation(int neighborhood_case, 
                        double llvel, 
                        double lvel, 
                        double cvel, 
                        double rvel, 
                        double rrvel, 
                        double &d2vdx2) 
{

    double v1, v2, v3, v4; // velocity definition see sketch

    if (neighborhood_case == 1) {
        // d2vdx2 = -8.0 * v2;
        d2vdx2 = -8.0 * cvel;
    }
    else if (neighborhood_case == 2) {
        // d2vdx2 = 8.0 / 3.0 * v3 - 16.0 / 3.0 * v2;
        d2vdx2 = 8.0 / 3.0 * rvel - 16.0 / 3.0 * cvel;
    }
    else if (neighborhood_case == 3) {
        // d2vdx2 = 8.0 / 3.0 * v2 - 16.0 / 3.0 * v3;
        d2vdx2 = 8.0 / 3.0 * lvel - 16.0 / 3.0 * cvel;
    }
    else if (neighborhood_case == 4) {
        // d2vdx2 = 2.0 * v3 - 5.0 * v2 - 1.0 / 5.0 * v4;
        d2vdx2 = 2.0 * rvel - 5.0 * cvel - 1.0 / 5.0 * rrvel;
    }
    else if (neighborhood_case == 5) {
        // d2vdx2 = 2.0 * v2 - 1.0 / 5.0 * v1 - 5.0 * v3;
        d2vdx2 = 2.0 * lvel - 1.0 / 5.0 * llvel - 5.0 * cvel;
    }
    else if (neighborhood_case == 6) {
        // d2vdx2 = v1 - 2.0 * v2 + v3;
        d2vdx2 = lvel - 2.0 * cvel + rvel;
    }
    return;
}


void split_number(  int number, 
                    int &case_x, 
                    int &case_y, 
                    int &case_z) 
{

    // three-digit number is split in three single numbers
    case_z = number%10;
    number /= 10;
    case_y = number%10;
    number /= 10;
    case_x = number%10;
    number /= 10;

    return;
}


void central_diff5( bool*** proc_geom, 
                    int* new_proc_size, 
                    double*** press, 
                    double*** vel_x, 
                    double*** vel_y, 
                    double*** vel_z, 
                    int*** voxel_neighborhood, 
                    unsigned int solver)
{
    int i, j, k;
    double dp_dx, dp_dy, dp_dz;             // difference quotient pressure
    double velx_sum = 0.0;                  // storing sum of all fluid vox velocites 
    double vely_sum = 0.0;                  // storing sum of all fluid vox velocites
    double velz_sum = 0.0;                  // storing sum of all fluid vox velocites
    double v_sum    = 0.0;                  // storing sum of all fluid vox velocites

    double d2vxdx2, d2vxdy2, d2vxdz2;       // Storing second order derivatives
    double d2vydx2, d2vydy2, d2vydz2;       // Storing second order derivatives
    double d2vzdx2, d2vzdy2, d2vzdz2;       // Storing second order derivatives

    int case_x, case_y, case_z;             // distingush neighborhood

    int lim = 2;                            // dependend on halo size

    for (i = lim; i < new_proc_size[2] - lim; i++){
        for (j = lim; j < new_proc_size[1] - lim; j++){
            for (k = lim; k < new_proc_size[0] - lim; k++){
                if (proc_geom[i][j][k] == 0){// only proceed if voxel {i, j, k} is element \omega^f
                    // cases for pressure gradient computation

                    if (proc_geom[i-1][j][k] == 0){
                        dp_dz = - press[i-1][j][k] + press[i][j][k];  
                    }
                    else {
                        dp_dz = 0.;
                    }

                    if (proc_geom[i][j-1][k] == 0){
                        dp_dy = - press[i][j-1][k] + press[i][j][k];  
                    }
                    else {
                        dp_dy = 0.;
                    }

                    if (proc_geom[i][j][k-1] == 0){
                        dp_dx = - press[i][j][k-1] + press[i][j][k]; 
                    }
                    else {
                        dp_dx = 0.;
                    }

                    // translate three-digit number in single cases
                    split_number(voxel_neighborhood[i][j][k], case_x, case_y, case_z);

                    d2vxdx2 = (vel_x[i][j][k-1] - 2*vel_x[i][j][k] + vel_x[i][j][k+1]); // flow and derivation in x-direction
                    d2vydy2 = (vel_y[i][j-1][k] - 2*vel_y[i][j][k] + vel_y[i][j+1][k]); // flow and derivation in y-direction
                    d2vzdz2 = (vel_z[i-1][j][k] - 2*vel_z[i][j][k] + vel_z[i+1][j][k]); // flow and derivation in z-direction

                    // perpendicular direction x
                    get_2nd_derivation(case_x, vel_y[i][j][k-2], vel_y[i][j][k-1], vel_y[i][j][k], vel_y[i][j][k+1], vel_y[i][j][k+2], d2vydx2);
                    get_2nd_derivation(case_x, vel_z[i][j][k-2], vel_z[i][j][k-1], vel_z[i][j][k], vel_z[i][j][k+1], vel_z[i][j][k+2], d2vzdx2);

                    // perpendicular direction y
                    get_2nd_derivation(case_y, vel_x[i][j-2][k], vel_x[i][j-1][k], vel_x[i][j][k], vel_x[i][j+1][k], vel_x[i][j+2][k], d2vxdy2);
                    get_2nd_derivation(case_y, vel_z[i][j-2][k], vel_z[i][j-1][k], vel_z[i][j][k], vel_z[i][j+1][k], vel_z[i][j+2][k], d2vzdy2);

                    // perpendicular direction y
                    get_2nd_derivation(case_z, vel_x[i-2][j][k], vel_x[i-1][j][k], vel_x[i][j][k], vel_x[i+1][j][k], vel_x[i+2][j][k], d2vxdz2);
                    get_2nd_derivation(case_z, vel_y[i-2][j][k], vel_y[i-1][j][k], vel_y[i][j][k], vel_y[i+1][j][k], vel_y[i+2][j][k], d2vydz2);

                    velx_sum = d2vxdx2 + d2vxdy2 + d2vxdz2;
                    vely_sum = d2vydy2 + d2vydz2 + d2vydx2;
                    velz_sum = d2vzdz2 + d2vzdx2 + d2vzdy2;
                    
                    // if solver is 1 (basically jacobi for the pressure) the pressure is computed with not updated vel_*[][][]
                    if (solver == 1){
                        v_sum = (vel_z[i+1][j][k] - vel_z[i][j][k] + vel_y[i][j+1][k] - vel_y[i][j][k] + vel_x[i][j][k+1] - vel_x[i][j][k]);
                    }
                    if (solver == 1 || solver == 2){
                        if (proc_geom[i][j][k-1] == 1){
                            vel_x[i][j][k] = 0.0;
                        }
                        else {
                            vel_x[i][j][k] = vel_x[i][j][k] + (-dp_dx + velx_sum / Re) * dt;
                        }
                        if (proc_geom[i][j-1][k] == 1){
                            vel_y[i][j][k] = 0.0;
                        }
                        else {
                            vel_y[i][j][k] = vel_y[i][j][k] + (-dp_dy + vely_sum / Re) * dt;
                        }
                        if (proc_geom[i-1][j][k] == 1){
                            vel_z[i][j][k] = 0.0;
                        }
                        else {
                            vel_z[i][j][k] = vel_z[i][j][k] + (-dp_dz + velz_sum / Re) * dt;
                        }
                    }
                    else if (solver == 3){
                        vel_x[i][j][k] = (1.0 - omega_sor) * vel_x[i][j][k] + omega_sor * (-dp_dx + velx_sum / Re) * dt;
                        vel_y[i][j][k] = (1.0 - omega_sor) * vel_y[i][j][k] + omega_sor * (-dp_dy + vely_sum / Re) * dt;
                        vel_z[i][j][k] = (1.0 - omega_sor) * vel_z[i][j][k] + omega_sor * (-dp_dz + velz_sum / Re) * dt;
                    }
                    // if solver is 2, 3 the pressure is computed with updated vel_*[][][]
                    if (solver == 2 || solver == 3){
                        v_sum = (vel_z[i+1][j][k] - vel_z[i][j][k] + vel_y[i][j+1][k] - vel_y[i][j][k] + vel_x[i][j][k+1] - vel_x[i][j][k]);
                    }
                    if (solver == 1 || solver == 2 || solver == 3){
                        press[i][j][k] = press[i][j][k] - c2 * v_sum * dt;
                    }
                    // else if (solver == 3){
                    //  press[i][j][k] = (1.0 - omega_sor) * press[i][j][k] - c2 * omega_sor * v_sum * dt;                      
                    // }
                }
            }
        }
    }
}