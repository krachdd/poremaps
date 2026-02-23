# SPDX-License-Identifier: MIT
#
# Copyright 2024 David Krach, Matthias Ruf
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

"""Shared pytest fixtures for the POREMAPS test suite."""

import numpy as np
import pytest

from helpers import make_proc_geom


def pytest_addoption(parser):
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Also run integration tests that invoke the POREMAPS binary.",
    )


@pytest.fixture
def run_integration(request):
    return request.config.getoption("--run-integration")


# ---------------------------------------------------------------------------
# Small geometry fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def full_fluid_geom():
    """8×8×8 all-fluid interior wrapped in solid halos; shape (12, 12, 12)."""
    interior = np.zeros((8, 8, 8), dtype=bool)
    return make_proc_geom(interior)


@pytest.fixture
def channel_z_geom():
    """10×6×6 channel open in z; solid walls at x=0/x_max and y=0/y_max.

    Interior shape (nz=10, ny=6, nx=6); fluid occupies y=1..4, x=1..4.
    Full array shape (14, 10, 10) after halos.
    """
    interior = np.zeros((10, 6, 6), dtype=bool)
    interior[:, 0,  :] = True   # solid y-wall (bottom)
    interior[:, -1, :] = True   # solid y-wall (top)
    interior[:, :,  0] = True   # solid x-wall (left)
    interior[:, :, -1] = True   # solid x-wall (right)
    return make_proc_geom(interior)


@pytest.fixture
def single_solid_voxel_geom():
    """6×6×6 mostly-fluid interior with one solid voxel at position (3,3,3).

    Full array shape (10, 10, 10) after halos.
    """
    interior = np.zeros((6, 6, 6), dtype=bool)
    interior[3, 3, 3] = True
    return make_proc_geom(interior)
