#!/usr/bin/python3

# Parallel Finite Difference Solver for Stokes Equations in Porous Media
# Copyright (C) 2023  David Krach, Matthias Ruf

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# author: dk, david.krach@mib.uni-stuttgart.de
# -- Header --------------------------------- #
import sys, os 
import numpy as np
import array
from pyevtk.hl import pointsToVTK as vtk
# -------------------------------------------- #

def import_rawfile(filename, size, dtype, order = 'F'):
    """
    filename : str
        name of the input file
    size : list [3]
        xsize, ysize, zsize (all three are int)
    dtype : data type of the input
        np.uint8 for input geometry
        np.uint32 for domain decomposition and voxel neigborhood
        np.float64 for pressure and velocity
        this may be platform dependent
    """
    # Total number of voxels
    count = size[0] * size[1] * size[2]
    # Open raw file
    rawfile = open(filename,'rb')
    # rawfile to numpy array
    array = np.fromfile(rawfile, dtype=dtype, count=count)
    # reshape array to 3d 
    array = np.reshape(array, (size[0], size[1], size[2]), order = order)
    rawfile.close()
    return array


# Sanity check
if len(sys.argv) != 6:
    raise ValueError('Usage: python3 field2vtu.py $FILENAME XSIZE YSIZE ZSIZE VOXELSIZE')
else:
    fname  = str(sys.argv[1])           # filename : str 
    xs     = int(sys.argv[2])           # xsize : int
    ys     = int(sys.argv[3])           # ysize : int 
    zs     = int(sys.argv[4])           # zsize : int 
    vs     = float(sys.argv[5])         # voxelsize : float (in meter)

x,y,z = np.mgrid[0:xs,0:ys,0:zs]*vs

# Shift positions by half voxelsize (to center of voxel)
x += 0.5*vs - 0.5*xs*vs
y += 0.5*vs - 0.5*ys*vs
z += 0.5*vs - 0.5*zs*vs

# Empty dict to save available data 
data = {}
# read geometry
fname_res = fname
if os.path.isfile(fname_res):
    domain = import_rawfile(fname_res, [xs, ys, zs], dtype = np.uint8)
    data['porespace'] = domain
else:
    raise ValueError('Geometry file not found.')

# read neighborhood
fname_res = f'voxel_neighborhood_{fname}'
if os.path.isfile(fname_res):
    voxnh = import_rawfile(fname_res, [xs, ys, zs], dtype = np.uint32)
    data['neighborhood'] = voxnh
    print(f'{fname_res} found.')
else:
    print(f'No voxel_neighborhood file found.')

# read neighborhood
fname_res = f'domain_decomp_{fname}'
if os.path.isfile(fname_res):
    decomp = import_rawfile(fname_res, [xs, ys, zs], dtype = np.uint32)
    data['decomposition'] = decomp
    print(f'{fname_res} found.')
else:
    print(f'No domain_decomp file found.')


# read press
fname_res = f'press_{fname}'
if os.path.isfile(fname_res):
    press = import_rawfile(fname_res, [xs, ys, zs], dtype = np.float64)
    data['pressure [-]'] = press
    print(f'{fname_res} found.')
else:
    print(f'No press file found.')

# read vels 
fname_res = f'velx_{fname}'
if os.path.isfile(fname_res):
    velx = import_rawfile(fname_res, [xs, ys, zs], dtype = np.float64)
    data['x_velocity [m/s]'] = velx
    print(f'{fname_res} found.')
else:
    print(f'No velx file found.')

fname_res = f'vely_{fname}'
if os.path.isfile(fname_res):
    vely = import_rawfile(fname_res, [xs, ys, zs], dtype = np.float64)
    data['y_velocity [m/s]'] = vely
    print(f'{fname_res} found.')
else:
    print(f'No vely file found.')

fname_res = f'velz_{fname}'
if os.path.isfile(fname_res):
    velz = import_rawfile(fname_res, [xs, ys, zs], dtype = np.float64)
    data['z_velocity [m/s]'] = velz
    print(f'{fname_res} found.')
else:
    print(f'No velz file found.')

# export results to one combined VTK File
filename_res = f'fields_{fname}'
vtk(filename_res.rsplit('.')[0], x, y, z, data = data)
