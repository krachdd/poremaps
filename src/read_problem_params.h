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
 * Purpose:    Header file to read_problem_params.cc
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/

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
                            int rank);


void read_parallel_params(  char* file_name, 
                            int* bc_method, 
                            int* dims,
                            int rank);