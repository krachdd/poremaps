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

"""Tests for geometry.cc functions: eval_geometry, get_dom_limits, get_proc_porosity."""

import numpy as np
import pytest

from helpers import (
    eval_geometry,
    get_dom_limits,
    get_proc_porosity,
    split_number,
    make_proc_geom,
    HALO,
)


class TestGetDomLimits:
    def test_single_rank_exact_division(self):
        size, dims, coords = [10, 8, 6], [1, 1, 1], [0, 0, 0]
        ps, nps, starts, ends = get_dom_limits(size, dims, coords)
        assert ps    == size
        assert nps   == [s + 4 for s in size]
        assert starts == [0, 0, 0]
        assert ends   == size

    def test_two_ranks_even_split(self):
        # 10 voxels into 2 ranks of 5 each; no remainder
        size, dims = [10, 10, 10], [2, 1, 1]
        ps0, nps0, st0, _ = get_dom_limits(size, dims, [0, 0, 0])
        ps1, nps1, st1, _ = get_dom_limits(size, dims, [1, 0, 0])
        assert ps0[0] == 5 and ps1[0] == 5
        assert st0[0] == 0 and st1[0] == 5
        assert nps0[0] == 9 and nps1[0] == 9

    def test_two_ranks_odd_split(self):
        # 11 voxels → rank 0 gets 6, rank 1 gets 5
        size, dims = [11, 8, 8], [2, 1, 1]
        ps0, _, st0, _ = get_dom_limits(size, dims, [0, 0, 0])
        ps1, _, st1, _ = get_dom_limits(size, dims, [1, 0, 0])
        assert ps0[0] + ps1[0] == 11
        assert ps0[0] == 6 and ps1[0] == 5
        # starts[1] accounts for rank 0's extra slice: 1*base + 1 = 6
        assert st1[0] == 6

    def test_starts_cover_full_domain_three_equal_ranks(self):
        # 12 voxels into 3 ranks of 4 each
        size, dims = [12, 12, 12], [3, 1, 1]
        starts = [get_dom_limits(size, dims, [c, 0, 0])[2][0] for c in range(3)]
        assert starts == [0, 4, 8]

    def test_new_proc_size_has_halos(self):
        size, dims = [10, 10, 10], [1, 1, 1]
        _, nps, _, _ = get_dom_limits(size, dims, [0, 0, 0])
        assert nps == [s + 4 for s in size]

    def test_three_uneven_ranks_sizes_sum_to_global(self):
        # 10 voxels into 3: base=3, remainder=1 → rank 0 gets 4, others get 3
        size, dims = [10, 8, 8], [3, 1, 1]
        sizes = [get_dom_limits(size, dims, [c, 0, 0])[0][0] for c in range(3)]
        assert sum(sizes) == size[0]
        assert sizes[0] == 4 and sizes[1] == 3 and sizes[2] == 3

    def test_y_and_z_also_split(self):
        # Decompose all three dimensions
        size, dims = [6, 9, 12], [2, 3, 4]
        total_x = sum(get_dom_limits(size, dims, [c, 0, 0])[0][0] for c in range(2))
        total_y = sum(get_dom_limits(size, dims, [0, c, 0])[0][1] for c in range(3))
        total_z = sum(get_dom_limits(size, dims, [0, 0, c])[0][2] for c in range(4))
        assert total_x == size[0]
        assert total_y == size[1]
        assert total_z == size[2]


class TestGetProcPorosity:
    def test_all_fluid(self, full_fluid_geom):
        assert get_proc_porosity(full_fluid_geom) == pytest.approx(1.0)

    def test_all_solid_interior(self):
        interior = np.ones((6, 6, 6), dtype=bool)
        geom = make_proc_geom(interior)
        assert get_proc_porosity(geom) == pytest.approx(0.0)

    def test_half_solid_half_fluid(self):
        interior = np.zeros((6, 6, 6), dtype=bool)
        interior[:, :, :3] = True   # solid on one x-half
        geom = make_proc_geom(interior)
        assert get_proc_porosity(geom) == pytest.approx(0.5)

    def test_halos_are_excluded(self):
        # Halos are solid by default; a fully-fluid interior still gives φ=1.
        interior = np.zeros((4, 4, 4), dtype=bool)
        geom = make_proc_geom(interior, halo_value=1)
        assert get_proc_porosity(geom) == pytest.approx(1.0)

    def test_channel_porosity(self, channel_z_geom):
        # Interior 10×6×6; walls at y=0, y=5, x=0, x=5 (4 fluid rows in each)
        # fluid fraction = (4*4) / (6*6) = 16/36
        assert get_proc_porosity(channel_z_geom) == pytest.approx(16 / 36)

    def test_single_solid_voxel(self, single_solid_voxel_geom):
        total  = 6 ** 3
        fluid  = total - 1
        assert get_proc_porosity(single_solid_voxel_geom) == pytest.approx(fluid / total)


class TestEvalGeometry:
    def test_solid_voxels_stay_minus_one(self):
        interior = np.zeros((6, 6, 6), dtype=bool)
        interior[3, 3, 3] = True   # one solid voxel
        geom  = make_proc_geom(interior)
        nhood = eval_geometry(geom)
        solid_mask = geom != 0
        assert (nhood[solid_mask] == -1).all()

    def test_fluid_voxels_have_nonnegative_code(self, full_fluid_geom):
        nhood    = eval_geometry(full_fluid_geom)
        interior = (slice(HALO, -HALO),) * 3
        fluid    = full_fluid_geom[interior] == 0
        assert (nhood[interior][fluid] >= 0).all()

    def test_isolated_fluid_voxel_all_directions_case1(self):
        # Single fluid voxel with solid on all sides → case 1 in every direction
        interior = np.ones((7, 7, 7), dtype=bool)
        interior[3, 3, 3] = False  # one fluid voxel
        geom  = make_proc_geom(interior)
        nhood = eval_geometry(geom)
        code  = nhood[3 + HALO, 3 + HALO, 3 + HALO]
        cx, cy, cz = split_number(code)
        assert cx == 1 and cy == 1 and cz == 1

    def test_fully_fluid_interior_voxel_all_directions_case6(self, full_fluid_geom):
        # A voxel well away from the halo boundary should have code 666
        nhood = eval_geometry(full_fluid_geom)
        nz, ny, nx = full_fluid_geom.shape
        mid = (nz // 2, ny // 2, nx // 2)
        assert nhood[mid] == 6 * 100 + 6 * 10 + 6

    def test_first_fluid_x_in_channel_is_near_wall(self, channel_z_geom):
        # First fluid voxel in x (interior x=1, array x=HALO+1=3) has:
        #   left  = solid wall (case_x has wall-left component)
        #   right = fluid, next-right = fluid  → case 4
        nhood = eval_geometry(channel_z_geom)
        nz, ny, _ = channel_z_geom.shape
        i = nz // 2
        j = ny // 2
        k = HALO + 1   # first fluid x-index in full array
        cx, _, _ = split_number(nhood[i, j, k])
        assert cx == 4

    def test_all_solid_geometry_all_minus_one(self):
        interior = np.ones((6, 6, 6), dtype=bool)
        geom  = make_proc_geom(interior)
        nhood = eval_geometry(geom)
        assert (nhood[HALO:-HALO, HALO:-HALO, HALO:-HALO] == -1).all()

    def test_two_fluid_voxels_adjacent_case2_and_case3(self):
        # Two adjacent fluid voxels in x surrounded by solid:
        #   left voxel: wall on left (lvox=1), fluid on right (rvox=0), wall beyond (rrvox=1) → case_x=2
        #   right voxel: fluid on left (lvox=0), wall on right (rvox=1), wall beyond (llvox=1) → case_x=3
        interior = np.ones((5, 5, 5), dtype=bool)
        interior[2, 2, 1] = False   # left fluid
        interior[2, 2, 2] = False   # right fluid
        geom  = make_proc_geom(interior)
        nhood = eval_geometry(geom)
        # left voxel: array index (2+HALO, 2+HALO, 1+HALO) = (4, 4, 3)
        cx_left,  _, _ = split_number(nhood[2 + HALO, 2 + HALO, 1 + HALO])
        # right voxel: array index (4, 4, 4)
        cx_right, _, _ = split_number(nhood[2 + HALO, 2 + HALO, 2 + HALO])
        assert cx_left  == 2
        assert cx_right == 3
