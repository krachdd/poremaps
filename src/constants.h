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
 * Purpose:    constants
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0

file_name_char      : Number of bytes reserved in the memory for saving a char series
pass                : Empty function to pass not yet specified functions
omega_sor           : Relaxation parameter for the SOR solver
Re                  : Reynolds number 
ndims               : Number of dimensions of Cartesian grid
c2                  : Dimensionless speed of sound
ndims               : Number of spatial dimensions
dt                  : Timestep size play with this number to decrease
                      convergence time
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

