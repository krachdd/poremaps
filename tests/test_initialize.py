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

"""Tests for initialize.cc functions: determine_comm_pattern, set_halos_initial."""

import numpy as np
import pytest

from helpers import determine_comm_pattern, set_halos_initial, HALO


class TestDetermineCommPattern:
    """Test bc[] array produced by determine_comm_pattern for each bc_method."""

    # ------------------------------------------------------------------
    # bc_method 0: fully periodic — all ranks, all directions → bc = 0
    # ------------------------------------------------------------------
    def test_method0_single_rank(self):
        assert determine_comm_pattern([1, 1, 1], [0, 0, 0], 0) == [0, 0, 0]

    def test_method0_interior_rank(self):
        assert determine_comm_pattern([4, 4, 4], [2, 2, 2], 0) == [0, 0, 0]

    def test_method0_boundary_rank_still_zero(self):
        assert determine_comm_pattern([4, 4, 4], [0, 0, 0], 0) == [0, 0, 0]

    # ------------------------------------------------------------------
    # bc_method 1: periodic z; slip walls in x and y
    # ------------------------------------------------------------------
    def test_method1_single_rank_in_all_dims(self):
        bc = determine_comm_pattern([1, 1, 1], [0, 0, 0], 1)
        assert bc[0] == 5 and bc[1] == 5 and bc[2] == 0

    def test_method1_z_is_always_zero(self):
        for coords in ([0, 0, 0], [1, 1, 1], [2, 2, 2]):
            bc = determine_comm_pattern([3, 3, 3], coords, 1)
            assert bc[2] == 0

    def test_method1_first_rank_in_x(self):
        bc = determine_comm_pattern([3, 1, 1], [0, 0, 0], 1)
        assert bc[0] == 1

    def test_method1_last_rank_in_x(self):
        bc = determine_comm_pattern([3, 1, 1], [2, 0, 0], 1)
        assert bc[0] == 3

    def test_method1_interior_rank_x_zero(self):
        bc = determine_comm_pattern([3, 1, 1], [1, 0, 0], 1)
        assert bc[0] == 0

    def test_method1_first_rank_in_y(self):
        bc = determine_comm_pattern([1, 4, 1], [0, 0, 0], 1)
        assert bc[1] == 1

    def test_method1_last_rank_in_y(self):
        bc = determine_comm_pattern([1, 4, 1], [0, 3, 0], 1)
        assert bc[1] == 3

    # ------------------------------------------------------------------
    # bc_method 2: periodic z; no-slip walls in x and y
    # ------------------------------------------------------------------
    def test_method2_single_rank(self):
        bc = determine_comm_pattern([1, 1, 1], [0, 0, 0], 2)
        assert bc[0] == 6 and bc[1] == 6 and bc[2] == 0

    def test_method2_first_rank_in_x(self):
        bc = determine_comm_pattern([2, 1, 1], [0, 0, 0], 2)
        assert bc[0] == 2

    def test_method2_last_rank_in_x(self):
        bc = determine_comm_pattern([2, 1, 1], [1, 0, 0], 2)
        assert bc[0] == 4

    def test_method2_interior_rank_zero(self):
        bc = determine_comm_pattern([3, 1, 1], [1, 0, 0], 2)
        assert bc[0] == 0

    def test_method2_z_is_always_zero(self):
        bc = determine_comm_pattern([1, 1, 3], [0, 0, 1], 2)
        assert bc[2] == 0

    # ------------------------------------------------------------------
    # bc_method 3: non-periodic; slip walls in all directions
    # ------------------------------------------------------------------
    def test_method3_single_rank_all(self):
        bc = determine_comm_pattern([1, 1, 1], [0, 0, 0], 3)
        assert bc == [5, 5, 5]

    def test_method3_first_rank_z(self):
        bc = determine_comm_pattern([1, 1, 3], [0, 0, 0], 3)
        assert bc[2] == 1

    def test_method3_last_rank_z(self):
        bc = determine_comm_pattern([1, 1, 3], [0, 0, 2], 3)
        assert bc[2] == 3

    def test_method3_interior_all_zero(self):
        bc = determine_comm_pattern([3, 3, 3], [1, 1, 1], 3)
        assert bc == [0, 0, 0]

    # ------------------------------------------------------------------
    # bc_method 4: non-periodic; no-slip x/y; no-slip z boundaries
    # ------------------------------------------------------------------
    def test_method4_single_rank_x_y_noslip(self):
        bc = determine_comm_pattern([1, 1, 1], [0, 0, 0], 4)
        assert bc[0] == 6 and bc[1] == 6

    def test_method4_single_rank_z_code(self):
        # Single rank in z gets code 5 (as per C++ code)
        bc = determine_comm_pattern([1, 1, 1], [0, 0, 0], 4)
        assert bc[2] == 5

    def test_method4_first_rank_z(self):
        bc = determine_comm_pattern([1, 1, 3], [0, 0, 0], 4)
        assert bc[2] == 2

    def test_method4_last_rank_z(self):
        bc = determine_comm_pattern([1, 1, 3], [0, 0, 2], 4)
        assert bc[2] == 4

    def test_method4_interior_all_zero(self):
        bc = determine_comm_pattern([3, 3, 3], [1, 1, 1], 4)
        assert bc == [0, 0, 0]

    # ------------------------------------------------------------------
    # All interior ranks are unaffected by bc_method (remain 0)
    # ------------------------------------------------------------------
    def test_interior_rank_unaffected_by_bc_method(self):
        for method in range(5):
            bc = determine_comm_pattern([5, 5, 5], [2, 2, 2], method)
            assert bc == [0, 0, 0], f"bc_method={method} changed interior rank bc"

    # ------------------------------------------------------------------
    # Default initialisation: bc starts at 0 for all ranks
    # ------------------------------------------------------------------
    def test_default_is_zero_before_overrides(self):
        # bc_method 0 leaves everything at 0
        bc = determine_comm_pattern([4, 4, 4], [1, 2, 2], 0)
        assert bc == [0, 0, 0]


class TestSetHalosInitial:
    def _ones(self, shape=(12, 12, 12)):
        return np.ones(shape, dtype=float)

    def test_zero_bc_leaves_array_unchanged(self):
        var = self._ones()
        ref = var.copy()
        set_halos_initial(var, [0, 0, 0])
        np.testing.assert_array_equal(var, ref)

    def test_slip_bcs_leave_array_unchanged(self):
        # bc values 1, 3, 5 are slip — they must NOT zero any halo layers
        for slip_bc in (
            [1, 0, 0], [3, 0, 0], [5, 0, 0],
            [0, 1, 0], [0, 3, 0], [0, 5, 0],
            [0, 0, 1], [0, 0, 3], [0, 0, 5],
            [5, 5, 5],
        ):
            var = self._ones()
            ref = var.copy()
            set_halos_initial(var, slip_bc)
            assert (var == ref).all(), f"bc={slip_bc} unexpectedly modified array"

    # X halos (axis 2, last dimension)
    def test_no_slip_x_left_zeros_columns_0_and_1(self):
        var = self._ones()
        set_halos_initial(var, [2, 0, 0])
        assert (var[:, :, 0] == 0.0).all()
        assert (var[:, :, 1] == 0.0).all()
        assert (var[:, :, -1] == 1.0).all()  # right side untouched

    def test_no_slip_x_right_zeros_last_two_columns(self):
        var = self._ones()
        set_halos_initial(var, [4, 0, 0])
        assert (var[:, :, -1] == 0.0).all()
        assert (var[:, :, -2] == 0.0).all()
        assert (var[:, :, 0] == 1.0).all()   # left side untouched

    def test_no_slip_x_both_zeros_all_x_halos(self):
        var = self._ones()
        set_halos_initial(var, [6, 0, 0])
        assert (var[:, :, 0] == 0.0).all()
        assert (var[:, :, 1] == 0.0).all()
        assert (var[:, :, -1] == 0.0).all()
        assert (var[:, :, -2] == 0.0).all()

    # Y halos (axis 1, middle dimension)
    def test_no_slip_y_left_zeros_rows_0_and_1(self):
        var = self._ones()
        set_halos_initial(var, [0, 2, 0])
        assert (var[:, 0, :] == 0.0).all()
        assert (var[:, 1, :] == 0.0).all()
        assert (var[:, -1, :] == 1.0).all()

    def test_no_slip_y_both_zeros_all_y_halos(self):
        var = self._ones()
        set_halos_initial(var, [0, 6, 0])
        assert (var[:, 0, :] == 0.0).all()
        assert (var[:, 1, :] == 0.0).all()
        assert (var[:, -1, :] == 0.0).all()
        assert (var[:, -2, :] == 0.0).all()

    # Z halos (axis 0, first dimension)
    def test_no_slip_z_left_zeros_planes_0_and_1(self):
        var = self._ones()
        set_halos_initial(var, [0, 0, 2])
        assert (var[0, :, :] == 0.0).all()
        assert (var[1, :, :] == 0.0).all()
        assert (var[-1, :, :] == 1.0).all()

    def test_no_slip_z_both_zeros_all_z_halos(self):
        var = self._ones()
        set_halos_initial(var, [0, 0, 6])
        assert (var[0, :, :] == 0.0).all()
        assert (var[1, :, :] == 0.0).all()
        assert (var[-1, :, :] == 0.0).all()
        assert (var[-2, :, :] == 0.0).all()

    def test_interior_never_zeroed_even_with_full_no_slip(self):
        var = self._ones()
        set_halos_initial(var, [6, 6, 6])
        interior = var[HALO:-HALO, HALO:-HALO, HALO:-HALO]
        assert (interior == 1.0).all()

    def test_all_no_slip_zeros_exactly_halo_layers(self):
        # With bc=[6,6,6], exactly the HALO=2 outer layers on each face should be 0.
        var = self._ones(shape=(10, 10, 10))
        set_halos_initial(var, [6, 6, 6])
        # Inner 6×6×6 block should be 1.0
        assert (var[2:-2, 2:-2, 2:-2] == 1.0).all()
        # Corner/edge/face halos should be 0.0 (spot checks)
        assert var[0, 0, 0] == 0.0
        assert var[1, 5, 5] == 0.0
        assert var[9, 9, 9] == 0.0
