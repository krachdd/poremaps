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

"""Tests for split_number and get_2nd_derivation (src/solver.cc)."""

import pytest

from helpers import split_number, get_2nd_derivation


class TestSplitNumber:
    def test_zero(self):
        assert split_number(0) == (0, 0, 0)

    def test_single_digit_z(self):
        assert split_number(6) == (0, 0, 6)

    def test_single_digit_y(self):
        assert split_number(30) == (0, 3, 0)

    def test_single_digit_x(self):
        assert split_number(200) == (2, 0, 0)

    def test_all_ones(self):
        assert split_number(111) == (1, 1, 1)

    def test_all_sixes(self):
        assert split_number(666) == (6, 6, 6)

    def test_mixed(self):
        assert split_number(246) == (2, 4, 6)
        assert split_number(135) == (1, 3, 5)

    def test_round_trip(self):
        """Encoding cases as x*100+y*10+z and decoding must recover originals."""
        for x in range(0, 7):
            for y in range(0, 7):
                for z in range(0, 7):
                    code = x * 100 + y * 10 + z
                    assert split_number(code) == (x, y, z)

    def test_only_case_z(self):
        cx, cy, cz = split_number(1)
        assert (cx, cy, cz) == (0, 0, 1)

    def test_only_case_x(self):
        cx, cy, cz = split_number(600)
        assert (cx, cy, cz) == (6, 0, 0)


class TestGet2ndDerivation:
    """Verify each neighbourhood case against the finite-difference formulae."""

    # ------------------------------------------------------------------
    # Case 1: both neighbours are walls  →  -8 * cvel
    # ------------------------------------------------------------------
    def test_case1_unit_centre(self):
        assert get_2nd_derivation(1, 0.0, 0.0, 1.0, 0.0, 0.0) == pytest.approx(-8.0)

    def test_case1_zero_centre(self):
        assert get_2nd_derivation(1, 5.0, 5.0, 0.0, 5.0, 5.0) == pytest.approx(0.0)

    def test_case1_ignores_neighbours(self):
        # llvel/lvel/rvel/rrvel are irrelevant for case 1
        r1 = get_2nd_derivation(1, 99.0, 99.0, 2.0, 99.0, 99.0)
        assert r1 == pytest.approx(-16.0)

    # ------------------------------------------------------------------
    # Case 2: wall on left, one fluid cell to the right
    # ------------------------------------------------------------------
    def test_case2_formula(self):
        result   = get_2nd_derivation(2, 0.0, 0.0, 1.0, 3.0, 0.0)
        expected = 8.0 / 3.0 * 3.0 - 16.0 / 3.0 * 1.0
        assert result == pytest.approx(expected)

    def test_case2_zero_result(self):
        # 8/3 * rvel == 16/3 * cvel  →  rvel = 2 * cvel
        assert get_2nd_derivation(2, 0.0, 0.0, 1.0, 2.0, 0.0) == pytest.approx(0.0)

    # ------------------------------------------------------------------
    # Case 3: wall on right, one fluid cell to the left
    # ------------------------------------------------------------------
    def test_case3_formula(self):
        result   = get_2nd_derivation(3, 0.0, 3.0, 1.0, 0.0, 0.0)
        expected = 8.0 / 3.0 * 3.0 - 16.0 / 3.0 * 1.0
        assert result == pytest.approx(expected)

    def test_case2_case3_mirror_symmetry(self):
        """Swapping left and right with equal values gives the same result."""
        r2 = get_2nd_derivation(2, 0.0, 0.0, 2.0, 4.0, 0.0)
        r3 = get_2nd_derivation(3, 0.0, 4.0, 2.0, 0.0, 0.0)
        assert r2 == pytest.approx(r3)

    # ------------------------------------------------------------------
    # Case 4: wall on left, two fluid cells to the right
    # ------------------------------------------------------------------
    def test_case4_formula(self):
        result   = get_2nd_derivation(4, 0.0, 0.0, 1.0, 2.0, 3.0)
        expected = 2.0 * 2.0 - 5.0 * 1.0 - 1.0 / 5.0 * 3.0
        assert result == pytest.approx(expected)

    # ------------------------------------------------------------------
    # Case 5: wall on right, two fluid cells to the left
    # ------------------------------------------------------------------
    def test_case5_formula(self):
        result   = get_2nd_derivation(5, 3.0, 2.0, 1.0, 0.0, 0.0)
        expected = 2.0 * 2.0 - 1.0 / 5.0 * 3.0 - 5.0 * 1.0
        assert result == pytest.approx(expected)

    def test_case4_case5_mirror_symmetry(self):
        """Mirroring stencil values between cases 4 and 5 gives equal results."""
        r4 = get_2nd_derivation(4, 0.0, 0.0, 2.0, 4.0, 6.0)
        r5 = get_2nd_derivation(5, 6.0, 4.0, 2.0, 0.0, 0.0)
        assert r4 == pytest.approx(r5)

    # ------------------------------------------------------------------
    # Case 6: fluid on both sides — standard 3-point central difference
    # ------------------------------------------------------------------
    def test_case6_quadratic_field(self):
        # v(x) = x²: v(-1)=1, v(0)=0, v(1)=1  →  d²v/dx² = 2
        assert get_2nd_derivation(6, 0.0, 1.0, 0.0, 1.0, 0.0) == pytest.approx(2.0)

    def test_case6_linear_field_zero_laplacian(self):
        # v(x) = x: v(-1)=-1, v(0)=0, v(1)=1  →  d²v/dx² = 0
        assert get_2nd_derivation(6, 0.0, -1.0, 0.0, 1.0, 0.0) == pytest.approx(0.0)

    def test_case6_constant_field_zero_laplacian(self):
        assert get_2nd_derivation(6, 5.0, 5.0, 5.0, 5.0, 5.0) == pytest.approx(0.0)

    def test_case6_ignores_second_neighbours(self):
        # llvel and rrvel are unused in case 6
        r1 = get_2nd_derivation(6, 99.0, 1.0, 2.0, 3.0, 99.0)
        r2 = get_2nd_derivation(6, -1.0, 1.0, 2.0, 3.0, -1.0)
        assert r1 == pytest.approx(r2)

    # ------------------------------------------------------------------
    # Error handling
    # ------------------------------------------------------------------
    def test_invalid_case_zero_raises(self):
        with pytest.raises(ValueError):
            get_2nd_derivation(0, 1.0, 1.0, 1.0, 1.0, 1.0)

    def test_invalid_case_seven_raises(self):
        with pytest.raises(ValueError):
            get_2nd_derivation(7, 1.0, 1.0, 1.0, 1.0, 1.0)
