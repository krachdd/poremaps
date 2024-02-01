#!/usr/bin/python3

# Copyright 2024 David Krach, Matthias Ruf

# Permission is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the “Software”), to deal in 
# the Software without restriction, including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS IN THE SOFTWARE.

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
