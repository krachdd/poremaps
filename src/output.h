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
 * Purpose:    Header for output.cc
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/

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
                        double voxelsize);


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
                        MPI_Datatype etype);


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
                                MPI_Datatype etype);



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
                    double tps);