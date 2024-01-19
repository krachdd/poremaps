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
 * Purpose:    Header for solver.cc
 *
 * Contents:   
 * 
 * Changelog:  Version: 1.0.0
 ***********************************************************************/


void split_number(  int number, 
                    int& case_x, 
                    int& case_y, 
                    int& case_z);


void get_2nd_derivation(int neighborhood_case, 
                        double llvel, 
                        double lvel, 
                        double cvel, 
                        double rvel, 
                        double rrvel, 
                        double& d2vdx2);


void central_diff5( bool*** proc_geom, 
                    int* new_proc_size, 
                    double*** press, 
                    double*** vel_x, 
                    double*** vel_y, 
                    double*** vel_z, 
                    int*** voxel_neighborhood, 
                    unsigned int solver);