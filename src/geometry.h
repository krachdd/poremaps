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
 * Purpose:    Evaluate Geometry. Neighbourhood. Read input geometry file (.raw)
 *
 * Contents:   Header for geometry.cc
 * 
 * Changelog:  Version: 1.0.0
 * 
 ***********************************************************************/

void read_geometry(     const char* GEOM_FILE, 
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
                        int* neighbors);

void get_dom_limits(    int* size, 
                        int* proc_size, 
                        int* new_proc_size, 
                        int* starts, 
                        int* ends, 
                        int* dims, 
                        int rank, 
                        int* cart_coords);

void eval_geometry(     bool*** proc_geom, 
                        int* new_proc_size, 
                        int*** voxel_neighborhood);

double get_proc_porosity(   bool*** proc_geom, 
                            int * new_proc_size);

double get_porosity(    bool*** proc_geom, 
                        int * new_proc_size, 
                        MPI_Comm comm_cart);
