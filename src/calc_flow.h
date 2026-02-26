/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright 2024-2026 David Krach, Matthias Ruf

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
 * Purpose:    Header for calc_flow.cc
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/

void calc_flow( bool*** proc_geom, 
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
                char* log_file);
