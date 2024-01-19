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
 * Purpose:    Header for evaluate.cc
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 * 
***********************************************************************/


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
                            double porosity);


double compute_convergence( bool*** proc_geom, 
                            int* size, 
                            int* new_proc_size, 
                            double*** vel_z, 
                            double* permeability,
                            MPI_Comm comm_cart);