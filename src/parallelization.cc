/************************************************************************

Parallel Finite Difference Solver for Stokes Equations in Porous Media
Copyright 2024 David Krach, Matthias Ruf

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
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
 * Purpose:    void communicate_geom_halos:
 *             used once, communicate halos geometry
 *
 *             void communicate_halos:
 *             pressure and velocity
 *             all 3 directions, if no communication is
 *             needed due to BCs or small number of ranks
 *             MPI_PROC_NULL means that there is no neighbour
 *             (defined as -2 in MPI4)
 * Contents:
 *
 * Changelog:  Version: 1.0.0

              krach 2021/01/14 added different cases for bc
                                added pindicator since pressure-bc have
                                to be treated different than velocity bc
                                in case of no slip conditions
                                if pindicator == 0 --> var is pressure
                                if pindicator == 1 --> var is velocity

              krach 2026: replaced element-by-element MPI_Sendrecv loops
                          with pack/send/unpack pattern.
                          Each direction now requires 2 MPI_Sendrecv calls
                          (6 total) instead of O(ny*nz) calls.
                          Also fixed MPI derived datatype leak in
                          communicate_geom_halos.
 ***********************************************************************/


// HEADER ***************************************************************
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "constants.h"

// **********************************************************************

void communicate_geom_halos(bool*** var,
                            int* new_proc_size,
                            int* neighbors,
                            MPI_Datatype etype,
                            MPI_Comm comm_cart)
{
    int i, j, k;
    MPI_Status status;

    int nx = new_proc_size[0];
    int ny = new_proc_size[1];
    int nz = new_proc_size[2];

    // --- X direction ---
    // Send x=2,3 to left neighbor (neighbors[0]); recv from right (neighbors[1]) into x=nx-2,nx-1
    // Send x=nx-4,nx-3 to right neighbor (neighbors[1]); recv from left (neighbors[0]) into x=0,1
    {
        int sz = nz * ny * 2;
        bool *sl = new bool[sz], *sr = new bool[sz];
        bool *rl = new bool[sz], *rr = new bool[sz];

        for (i = 0; i < nz; i++)
            for (j = 0; j < ny; j++) {
                int idx = 2*(i*ny + j);
                sl[idx]   = var[i][j][2];     sl[idx+1] = var[i][j][3];
                sr[idx]   = var[i][j][nx-4];  sr[idx+1] = var[i][j][nx-3];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[0], 0,
                     rr, sz, etype, neighbors[1], 0, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[1], 1,
                     rl, sz, etype, neighbors[0], 1, comm_cart, &status);

        if (neighbors[1] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (j = 0; j < ny; j++) {
                    int idx = 2*(i*ny + j);
                    var[i][j][nx-2] = rr[idx]; var[i][j][nx-1] = rr[idx+1];
                }

        if (neighbors[0] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (j = 0; j < ny; j++) {
                    int idx = 2*(i*ny + j);
                    var[i][j][0] = rl[idx]; var[i][j][1] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }

    // --- Y direction ---
    // Send y=2,3 to bottom neighbor (neighbors[2]); recv from top (neighbors[3]) into y=ny-2,ny-1
    // Send y=ny-4,ny-3 to top neighbor (neighbors[3]); recv from bottom (neighbors[2]) into y=0,1
    {
        int sz = nz * nx * 2;
        bool *sl = new bool[sz], *sr = new bool[sz];
        bool *rl = new bool[sz], *rr = new bool[sz];

        for (i = 0; i < nz; i++)
            for (k = 0; k < nx; k++) {
                int idx = 2*(i*nx + k);
                sl[idx]   = var[i][2][k];     sl[idx+1] = var[i][3][k];
                sr[idx]   = var[i][ny-4][k];  sr[idx+1] = var[i][ny-3][k];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[2], 2,
                     rr, sz, etype, neighbors[3], 2, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[3], 3,
                     rl, sz, etype, neighbors[2], 3, comm_cart, &status);

        if (neighbors[3] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(i*nx + k);
                    var[i][ny-2][k] = rr[idx]; var[i][ny-1][k] = rr[idx+1];
                }

        if (neighbors[2] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(i*nx + k);
                    var[i][0][k] = rl[idx]; var[i][1][k] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }

    // --- Z direction ---
    // Send z=2,3 to front neighbor (neighbors[4]); recv from back (neighbors[5]) into z=nz-2,nz-1
    // Send z=nz-4,nz-3 to back neighbor (neighbors[5]); recv from front (neighbors[4]) into z=0,1
    {
        int sz = ny * nx * 2;
        bool *sl = new bool[sz], *sr = new bool[sz];
        bool *rl = new bool[sz], *rr = new bool[sz];

        for (j = 0; j < ny; j++)
            for (k = 0; k < nx; k++) {
                int idx = 2*(j*nx + k);
                sl[idx]   = var[2][j][k];     sl[idx+1] = var[3][j][k];
                sr[idx]   = var[nz-4][j][k];  sr[idx+1] = var[nz-3][j][k];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[4], 4,
                     rr, sz, etype, neighbors[5], 4, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[5], 5,
                     rl, sz, etype, neighbors[4], 5, comm_cart, &status);

        if (neighbors[5] != MPI_PROC_NULL)
            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(j*nx + k);
                    var[nz-2][j][k] = rr[idx]; var[nz-1][j][k] = rr[idx+1];
                }

        if (neighbors[4] != MPI_PROC_NULL)
            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(j*nx + k);
                    var[0][j][k] = rl[idx]; var[1][j][k] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }
}


void communicate_halos( double*** var,
                        int* new_proc_size,
                        int* neighbors,
                        MPI_Datatype etype,
                        MPI_Comm comm_cart,
                        unsigned int* bc,
                        int pindicator)
{
    int i, j, k;
    MPI_Status status;

    int nx = new_proc_size[0];
    int ny = new_proc_size[1];
    int nz = new_proc_size[2];

    // --- X direction ---
    // Send x=2,3 to left neighbor (neighbors[0]); recv from right (neighbors[1]) into x=nx-2,nx-1
    // Send x=nx-4,nx-3 to right neighbor (neighbors[1]); recv from left (neighbors[0]) into x=0,1
    {
        int sz = nz * ny * 2;
        double *sl = new double[sz], *sr = new double[sz];
        double *rl = new double[sz], *rr = new double[sz];

        for (i = 0; i < nz; i++)
            for (j = 0; j < ny; j++) {
                int idx = 2*(i*ny + j);
                sl[idx]   = var[i][j][2];     sl[idx+1] = var[i][j][3];
                sr[idx]   = var[i][j][nx-4];  sr[idx+1] = var[i][j][nx-3];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[0], 0,
                     rr, sz, etype, neighbors[1], 0, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[1], 1,
                     rl, sz, etype, neighbors[0], 1, comm_cart, &status);

        if (neighbors[1] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (j = 0; j < ny; j++) {
                    int idx = 2*(i*ny + j);
                    var[i][j][nx-2] = rr[idx]; var[i][j][nx-1] = rr[idx+1];
                }

        if (neighbors[0] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (j = 0; j < ny; j++) {
                    int idx = 2*(i*ny + j);
                    var[i][j][0] = rl[idx]; var[i][j][1] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }

    // --- Y direction ---
    // Send y=2,3 to bottom neighbor (neighbors[2]); recv from top (neighbors[3]) into y=ny-2,ny-1
    // Send y=ny-4,ny-3 to top neighbor (neighbors[3]); recv from bottom (neighbors[2]) into y=0,1
    {
        int sz = nz * nx * 2;
        double *sl = new double[sz], *sr = new double[sz];
        double *rl = new double[sz], *rr = new double[sz];

        for (i = 0; i < nz; i++)
            for (k = 0; k < nx; k++) {
                int idx = 2*(i*nx + k);
                sl[idx]   = var[i][2][k];     sl[idx+1] = var[i][3][k];
                sr[idx]   = var[i][ny-4][k];  sr[idx+1] = var[i][ny-3][k];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[2], 2,
                     rr, sz, etype, neighbors[3], 2, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[3], 3,
                     rl, sz, etype, neighbors[2], 3, comm_cart, &status);

        if (neighbors[3] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(i*nx + k);
                    var[i][ny-2][k] = rr[idx]; var[i][ny-1][k] = rr[idx+1];
                }

        if (neighbors[2] != MPI_PROC_NULL)
            for (i = 0; i < nz; i++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(i*nx + k);
                    var[i][0][k] = rl[idx]; var[i][1][k] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }

    // --- Z direction ---
    // Send z=2,3 to front neighbor (neighbors[4]); recv from back (neighbors[5]) into z=nz-2,nz-1
    // Send z=nz-4,nz-3 to back neighbor (neighbors[5]); recv from front (neighbors[4]) into z=0,1
    {
        int sz = ny * nx * 2;
        double *sl = new double[sz], *sr = new double[sz];
        double *rl = new double[sz], *rr = new double[sz];

        for (j = 0; j < ny; j++)
            for (k = 0; k < nx; k++) {
                int idx = 2*(j*nx + k);
                sl[idx]   = var[2][j][k];     sl[idx+1] = var[3][j][k];
                sr[idx]   = var[nz-4][j][k];  sr[idx+1] = var[nz-3][j][k];
            }

        MPI_Sendrecv(sl, sz, etype, neighbors[4], 4,
                     rr, sz, etype, neighbors[5], 4, comm_cart, &status);
        MPI_Sendrecv(sr, sz, etype, neighbors[5], 5,
                     rl, sz, etype, neighbors[4], 5, comm_cart, &status);

        if (neighbors[5] != MPI_PROC_NULL)
            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(j*nx + k);
                    var[nz-2][j][k] = rr[idx]; var[nz-1][j][k] = rr[idx+1];
                }

        if (neighbors[4] != MPI_PROC_NULL)
            for (j = 0; j < ny; j++)
                for (k = 0; k < nx; k++) {
                    int idx = 2*(j*nx + k);
                    var[0][j][k] = rl[idx]; var[1][j][k] = rl[idx+1];
                }

        delete[] sl; delete[] sr; delete[] rl; delete[] rr;
    }
}
