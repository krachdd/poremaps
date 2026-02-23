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
 * Purpose:    constants
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 *
 * file_name_char  : max bytes for a filename string
 * pass            : no-op placeholder macro
 * omega_sor       : relaxation factor for the SOR solver (1.0 = Gauss-Seidel)
 * Re              : dimensionless viscosity (Reynolds-like parameter)
 * ndims           : number of spatial dimensions of the Cartesian grid
 * c2              : dimensionless speed-of-sound squared (artificial compressibility)
 * dt              : pseudo-time step size; reduce to improve stability
 ***********************************************************************/

#ifndef file_name_char
    #define file_name_char 400
#endif
#ifndef pass
    #define pass (void)0 
#endif
#ifndef omega_sor
    #define omega_sor 1.2 
#endif
#ifndef Re 
    #define Re 0.01
#endif
#ifndef c2 
    #define c2 1.5e6
#endif
#ifndef ndims
    #define ndims 3
#endif
#ifndef dt 
    #define dt 1.0e-4
#endif

