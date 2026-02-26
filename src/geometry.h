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
