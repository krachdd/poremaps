# SPDX-License-Identifier: MIT
#
# Copyright 2024-2026 David Krach, Matthias Ruf
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Pure-Python reference implementations of POREMAPS C++ kernels.

Used by the pytest suite to verify numerical logic independently of the
MPI-parallel C++ code.

Array conventions
-----------------
All 3-D arrays use numpy shape (nz+4, ny+4, nx+4), matching the C++
convention where the first index is z, second is y, third is x.  The
physical (inner) domain occupies ``[HALO:-HALO]`` in every dimension.
Geometry encoding: 0 = fluid, non-zero = solid.
"""

import numpy as np

# Physical constants — must match src/constants.h
Re       = 0.01
c2       = 1.5e6
dt       = 1.0e-4
omega_sor = 1.2
HALO     = 2   # ghost-layer thickness on each face


# ---------------------------------------------------------------------------
# solver.cc
# ---------------------------------------------------------------------------

def split_number(number):
    """Decode a 3-digit base-10 neighbourhood code into (case_x, case_y, case_z).

    Mirrors solver.cc::split_number().
    """
    case_z = number % 10
    number //= 10
    case_y = number % 10
    number //= 10
    case_x = number % 10
    return case_x, case_y, case_z


def get_2nd_derivation(case, llvel, lvel, cvel, rvel, rrvel):
    """Return the FD approximation to d²v/dx² for the given neighbourhood case.

    Mirrors solver.cc::get_2nd_derivation().

    Parameters
    ----------
    case : int
        Neighbourhood case 1–6.
    llvel, lvel, cvel, rvel, rrvel : float
        Stencil values: second-left, left, centre, right, second-right.
    """
    if case == 1:
        return -8.0 * cvel
    if case == 2:
        return 8.0 / 3.0 * rvel - 16.0 / 3.0 * cvel
    if case == 3:
        return 8.0 / 3.0 * lvel - 16.0 / 3.0 * cvel
    if case == 4:
        return 2.0 * rvel - 5.0 * cvel - 1.0 / 5.0 * rrvel
    if case == 5:
        return 2.0 * lvel - 1.0 / 5.0 * llvel - 5.0 * cvel
    if case == 6:
        return lvel - 2.0 * cvel + rvel
    raise ValueError(f"neighbourhood case must be 1–6, got {case}")


# ---------------------------------------------------------------------------
# geometry.cc
# ---------------------------------------------------------------------------

def get_dom_limits(size, dims, cart_coords):
    """Compute per-rank domain limits.

    Mirrors geometry.cc::get_dom_limits().

    Parameters
    ----------
    size        : sequence[3]  Global domain size [nx, ny, nz].
    dims        : sequence[3]  Cartesian decomposition [dx, dy, dz].
    cart_coords : sequence[3]  This rank's coordinates in the Cartesian grid.

    Returns
    -------
    proc_size     : list[3]  Owned voxels per direction (after extra-slice assignment).
    new_proc_size : list[3]  proc_size + 2*HALO.
    starts        : list[3]  Global start indices (based on base proc_size).
    ends          : list[3]  starts + base proc_size.
    """
    base     = [size[i] // dims[i] for i in range(3)]
    addslice = [size[i] - dims[i] * base[i] for i in range(3)]

    starts = [0, 0, 0]
    ends   = [0, 0, 0]
    for i in range(3):
        if cart_coords[i] == 0:
            cum = 0
        elif 0 < cart_coords[i] <= addslice[i]:
            cum = cart_coords[i]
        else:
            cum = addslice[i]
        starts[i] = cart_coords[i] * base[i] + cum
        ends[i]   = starts[i] + base[i]

    # Extra slice for early ranks (matches C++ order: ends computed before update)
    proc_size = [base[i] + (1 if cart_coords[i] < addslice[i] else 0)
                 for i in range(3)]
    new_proc_size = [p + 2 * HALO for p in proc_size]
    return proc_size, new_proc_size, starts, ends


def eval_geometry(proc_geom):
    """Encode voxel neighbourhood into a 3-digit base-10 code for each fluid voxel.

    Mirrors geometry.cc::eval_geometry() (Bentz et al. 2007).

    Parameters
    ----------
    proc_geom : ndarray, bool, shape (nz, ny, nx)
        0 = fluid, non-zero = solid (including halos).

    Returns
    -------
    nhood : ndarray, int32, same shape
        Neighbourhood code: 0 for unconstrained fluid, positive 3-digit code for
        near-wall fluid, -1 for solid voxels.  Only interior voxels are processed.
    """
    nhood = np.where(proc_geom == 0, 0, -1).astype(np.int32)
    nz, ny, nx = proc_geom.shape
    lim = HALO

    # (shift vector in [z,y,x], shift_factor)
    directions = [
        ((0, 0, 1), 100),   # x direction
        ((0, 1, 0), 10),    # y direction
        ((1, 0, 0), 1),     # z direction
    ]

    for (di, dj, dk), sf in directions:
        for i in range(lim, nz - lim):
            for j in range(lim, ny - lim):
                for k in range(lim, nx - lim):
                    if proc_geom[i, j, k] != 0:
                        continue  # skip solid voxels

                    lv  = 0 if proc_geom[i - di, j - dj, k - dk] == 0 else 1
                    rv  = 0 if proc_geom[i + di, j + dj, k + dk] == 0 else 1
                    llv = -1
                    rrv = -1
                    if lv == 0:
                        llv = 0 if proc_geom[i - 2*di, j - 2*dj, k - 2*dk] == 0 else 1
                    if rv == 0:
                        rrv = 0 if proc_geom[i + 2*di, j + 2*dj, k + 2*dk] == 0 else 1

                    if lv == 1 and rv == 1:
                        case = 1
                    elif lv == 1 and rv == 0 and rrv == 1:
                        case = 2
                    elif llv == 1 and lv == 0 and rv == 1:
                        case = 3
                    elif lv == 1 and rv == 0 and rrv == 0:
                        case = 4
                    elif llv == 0 and lv == 0 and rv == 1:
                        case = 5
                    elif lv == 0 and rv == 0:
                        case = 6
                    else:
                        case = 0  # unclassified; should not occur in valid geometry

                    nhood[i, j, k] += case * sf

    return nhood


def get_proc_porosity(proc_geom):
    """Return the porosity of the interior domain (halos excluded).

    Mirrors geometry.cc::get_proc_porosity().
    """
    interior   = proc_geom[HALO:-HALO, HALO:-HALO, HALO:-HALO]
    fluid_count = int(np.count_nonzero(interior == 0))
    vox_count   = interior.size
    return fluid_count / vox_count


# ---------------------------------------------------------------------------
# initialize.cc
# ---------------------------------------------------------------------------

def determine_comm_pattern(dims, cart_coords, bc_method):
    """Return bc[3] for the given rank in the Cartesian topology.

    Mirrors initialize.cc::determine_comm_pattern().

    bc values per direction:
      0: communicate both ways (interior rank)
      1: communicate +1, slip at -1
      2: communicate +1, no-slip at -1
      3: communicate -1, slip at +1
      4: communicate -1, no-slip at +1
      5: slip at both -1 and +1
      6: no-slip at both -1 and +1

    Parameters
    ----------
    dims        : sequence[3]  Cartesian topology dimensions [dx, dy, dz].
    cart_coords : sequence[3]  This rank's Cartesian coordinates.
    bc_method   : int          0–4 physical setup selector.

    Returns
    -------
    bc : list[3]
    """
    bc = [0, 0, 0]

    if bc_method == 0:
        pass  # fully periodic — bc stays 0

    elif bc_method in (1, 2):
        # Periodic z; slip (method 1) or no-slip (method 2) walls in x and y.
        lo   = 1 if bc_method == 1 else 2   # first-rank code
        hi   = 3 if bc_method == 1 else 4   # last-rank code
        both = 5 if bc_method == 1 else 6   # single-rank code
        bc[2] = 0  # z is always periodic
        for ax in (0, 1):
            if dims[ax] == 1:
                bc[ax] = both
            elif dims[ax] > 1:
                if cart_coords[ax] == 0:
                    bc[ax] = lo
                if cart_coords[ax] + 1 == dims[ax]:
                    bc[ax] = hi

    elif bc_method == 3:
        # Non-periodic; slip walls in all three directions.
        for ax in range(3):
            if dims[ax] == 1:
                bc[ax] = 5
            elif dims[ax] > 1:
                if cart_coords[ax] == 0:
                    bc[ax] = 1
                if cart_coords[ax] + 1 == dims[ax]:
                    bc[ax] = 3

    elif bc_method == 4:
        # Non-periodic; no-slip x/y walls; no-slip z boundaries.
        for ax in (0, 1):
            if dims[ax] == 1:
                bc[ax] = 6
            elif dims[ax] > 1:
                if cart_coords[ax] == 0:
                    bc[ax] = 2
                if cart_coords[ax] + 1 == dims[ax]:
                    bc[ax] = 4
        if dims[2] == 1:
            bc[2] = 5
        elif dims[2] > 1:
            if cart_coords[2] == 0:
                bc[2] = 2
            if cart_coords[2] + 1 == dims[2]:
                bc[2] = 4

    return bc


def set_halos_initial(var, bc):
    """Zero velocity halo layers that correspond to no-slip wall boundaries (in-place).

    Mirrors initialize.cc::set_halos_initial().

    Parameters
    ----------
    var : ndarray, float, shape (nz, ny, nx)
    bc  : sequence[3]  Boundary condition codes [x, y, z].
    """
    # X halos (last axis)
    if bc[0] in (2, 6):
        var[:, :, 0] = 0.0
        var[:, :, 1] = 0.0
    if bc[0] in (4, 6):
        var[:, :, -1] = 0.0
        var[:, :, -2] = 0.0

    # Y halos (middle axis)
    if bc[1] in (2, 6):
        var[:, 0, :] = 0.0
        var[:, 1, :] = 0.0
    if bc[1] in (4, 6):
        var[:, -1, :] = 0.0
        var[:, -2, :] = 0.0

    # Z halos (first axis)
    if bc[2] in (2, 6):
        var[0, :, :] = 0.0
        var[1, :, :] = 0.0
    if bc[2] in (4, 6):
        var[-1, :, :] = 0.0
        var[-2, :, :] = 0.0


# ---------------------------------------------------------------------------
# evaluation.cc
# ---------------------------------------------------------------------------

def compute_permeability_whole_domain(vel_x, vel_y, vel_z, proc_geom,
                                      voxelsize, porosity):
    """Compute whole-domain permeability components (single rank, no MPI).

    Mirrors the whole-domain section of evaluation.cc::compute_permeability().

    Returns
    -------
    wk13, wk23, wk33 : float  Permeability tensor components [m²].
    wmean_velz        : float  Mean z-velocity [m/s].
    wmax_velz         : float  Maximum z-velocity [m/s].
    """
    sl   = (slice(HALO, -HALO),) * 3
    fluid = proc_geom[sl] == 0

    if not fluid.any():
        return 0.0, 0.0, 0.0, 0.0, 0.0

    wflcount  = int(fluid.sum())
    wvelx_sum = float(vel_x[sl][fluid].sum())
    wvely_sum = float(vel_y[sl][fluid].sum())
    wvelz_sum = float(vel_z[sl][fluid].sum())
    wvelz_max = float(vel_z[sl][fluid].max())

    wk13 = (wvelx_sum / wflcount) * (1.0 / Re) * porosity * voxelsize ** 2
    wk23 = (wvely_sum / wflcount) * (1.0 / Re) * porosity * voxelsize ** 2
    wk33 = (wvelz_sum / wflcount) * (1.0 / Re) * porosity * voxelsize ** 2

    wmean_velz = wvelz_sum / wflcount * voxelsize
    wmax_velz  = wvelz_max * voxelsize

    return wk13, wk23, wk33, wmean_velz, wmax_velz


def compute_convergence(vel_z, proc_geom, permeability):
    """Compute the convergence criterion (single rank, no MPI).

    Mirrors evaluation.cc::compute_convergence().

    Parameters
    ----------
    vel_z        : ndarray, float, shape (nz, ny, nx)
    proc_geom    : ndarray, bool,  shape (nz, ny, nx)
    permeability : float  Permeability from the previous iteration.

    Returns
    -------
    conv     : float  Relative change |k_new - k_old| / k_new.
    new_perm : float  Updated permeability value.
    """
    sl    = (slice(HALO, -HALO),) * 3
    fluid = proc_geom[sl] == 0

    if not fluid.any():
        return 0.0, permeability

    velz_sum = float(vel_z[sl][fluid].sum())
    flcount  = int(fluid.sum())

    mean_velz            = velz_sum / flcount
    current_permeability = mean_velz / Re
    if current_permeability != 0.0:
        conv = abs(current_permeability - permeability) / current_permeability
    else:
        conv = 0.0

    return conv, current_permeability


# ---------------------------------------------------------------------------
# Geometry construction helpers
# ---------------------------------------------------------------------------

def make_proc_geom(geom_interior, halo_value=1):
    """Wrap an interior geometry array with HALO solid layers on every face.

    Parameters
    ----------
    geom_interior : ndarray, shape (nz, ny, nx)  0=fluid, non-zero=solid.
    halo_value    : scalar  Value to fill halo cells with (default 1 = solid).

    Returns
    -------
    ndarray, same dtype, shape (nz+4, ny+4, nx+4).
    """
    nz, ny, nx = geom_interior.shape
    out = np.full((nz + 2*HALO, ny + 2*HALO, nx + 2*HALO),
                  halo_value, dtype=geom_interior.dtype)
    out[HALO:-HALO, HALO:-HALO, HALO:-HALO] = geom_interior
    return out


def write_geometry(geom_interior, path):
    """Write a geometry to a binary .raw file (uint8, x-fastest / Fortran order).

    Parameters
    ----------
    geom_interior : ndarray, shape (nz, ny, nx)  0=fluid, 1=solid.
    path          : str or Path-like.
    """
    geom_interior.astype(np.uint8).ravel(order='C').tofile(path)
