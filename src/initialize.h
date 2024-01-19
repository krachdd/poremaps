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
 * Date:       2021/02
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    
 *
 * Contents:   Header for initialize.cc
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/ 

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
                            unsigned int* bc);

void determine_comm_pattern(int* dims, 
                            int rank, 
                            int* cart_coords, 
                            unsigned int* bc, 
                            int bc_method, 
                            int* neighbors);

void set_halos_initial( double*** var, 
                        int* new_proc_size, 
                        unsigned int* bc);
