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