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
 * Date:       2021/02
 * Contact:    david.krach@mib.uni-stuttgart.de
 *
 * Purpose:    Header for read_problem_params.cc
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