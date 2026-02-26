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

"""Tests for evaluation.cc functions: compute_permeability, compute_convergence."""

import numpy as np
import pytest

from helpers import (
    compute_permeability_whole_domain,
    compute_convergence,
    make_proc_geom,
    HALO,
    Re,
)


def _fluid_geom(nz=6, ny=6, nx=6):
    return make_proc_geom(np.zeros((nz, ny, nx), dtype=bool))


def _vel(value, shape):
    return np.full(shape, value, dtype=float)


class TestComputePermeabilityWholeDomain:
    def test_uniform_velz_gives_correct_k33(self):
        geom = _fluid_geom()
        s    = geom.shape
        dx   = 1e-6
        phi  = 1.0
        wk13, wk23, wk33, _, _ = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), _vel(1.0, s), geom, dx, phi
        )
        assert wk33 == pytest.approx((1.0 / Re) * phi * dx**2)
        assert wk13 == pytest.approx(0.0)
        assert wk23 == pytest.approx(0.0)

    def test_uniform_velx_gives_correct_k13(self):
        geom = _fluid_geom()
        s    = geom.shape
        dx   = 1e-6
        wk13, wk23, wk33, _, _ = compute_permeability_whole_domain(
            _vel(1.0, s), np.zeros(s), np.zeros(s), geom, dx, 1.0
        )
        assert wk13 == pytest.approx((1.0 / Re) * dx**2)
        assert wk23 == pytest.approx(0.0)
        assert wk33 == pytest.approx(0.0)

    def test_zero_velocity_gives_zero_permeability(self):
        geom = _fluid_geom()
        s    = geom.shape
        z    = np.zeros(s)
        wk13, wk23, wk33, wmean, wmax = compute_permeability_whole_domain(
            z, z, z, geom, 1e-6, 1.0
        )
        assert wk13 == pytest.approx(0.0)
        assert wk23 == pytest.approx(0.0)
        assert wk33 == pytest.approx(0.0)
        assert wmean == pytest.approx(0.0)
        assert wmax  == pytest.approx(0.0)

    def test_porosity_scales_permeability_linearly(self):
        geom  = _fluid_geom()
        s     = geom.shape
        vz    = _vel(1.0, s)
        dx    = 1e-6
        _, _, k33_p05, _, _ = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, dx, 0.5
        )
        _, _, k33_p10, _, _ = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, dx, 1.0
        )
        assert k33_p05 == pytest.approx(k33_p10 * 0.5)

    def test_voxelsize_scales_permeability_as_squared(self):
        geom = _fluid_geom()
        s    = geom.shape
        vz   = _vel(1.0, s)
        _, _, k33_1, _, _ = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, 1.0, 1.0
        )
        _, _, k33_2, _, _ = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, 2.0, 1.0
        )
        assert k33_2 == pytest.approx(k33_1 * 4.0)

    def test_mean_and_max_velz_scaled_by_voxelsize(self):
        geom  = _fluid_geom(4, 4, 4)
        s     = geom.shape
        vz    = np.zeros(s)
        vz[HALO:-HALO, HALO:-HALO, HALO:-HALO] = 3.0
        dx    = 2.0
        _, _, _, wmean, wmax = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, dx, 1.0
        )
        assert wmean == pytest.approx(3.0 * dx)
        assert wmax  == pytest.approx(3.0 * dx)

    def test_max_velz_is_maximum_not_mean(self):
        # Half of interior has vel=1, other half has vel=2 → max=2, mean=1.5
        geom  = _fluid_geom(4, 4, 4)
        s     = geom.shape
        vz    = np.zeros(s)
        nz, ny, nx = 4, 4, 4
        vz[HALO:HALO+nz//2, HALO:-HALO, HALO:-HALO] = 1.0
        vz[HALO+nz//2:-HALO, HALO:-HALO, HALO:-HALO] = 2.0
        _, _, _, wmean, wmax = compute_permeability_whole_domain(
            np.zeros(s), np.zeros(s), vz, geom, 1.0, 1.0
        )
        assert wmax  == pytest.approx(2.0)
        assert wmean == pytest.approx(1.5)

    def test_all_solid_returns_zeros(self):
        geom  = make_proc_geom(np.ones((6, 6, 6), dtype=bool))
        s     = geom.shape
        ones  = np.ones(s)
        wk13, wk23, wk33, wmean, wmax = compute_permeability_whole_domain(
            ones, ones, ones, geom, 1e-6, 0.0
        )
        assert wk13 == 0.0 and wk23 == 0.0 and wk33 == 0.0
        assert wmean == 0.0 and wmax == 0.0


class TestComputeConvergence:
    def test_zero_change_gives_zero_convergence(self):
        geom  = _fluid_geom()
        s     = geom.shape
        vel_z = _vel(1.0, s)
        # Set old permeability to exactly what the current field will give.
        old_perm = 1.0 / Re
        conv, new_perm = compute_convergence(vel_z, geom, old_perm)
        assert conv     == pytest.approx(0.0)
        assert new_perm == pytest.approx(old_perm)

    def test_permeability_updated_on_each_call(self):
        geom = _fluid_geom()
        s    = geom.shape
        _, p1 = compute_convergence(_vel(1.0, s), geom, 0.0)
        _, p2 = compute_convergence(_vel(2.0, s), geom, p1)
        assert p1 == pytest.approx(1.0 / Re)
        assert p2 == pytest.approx(2.0 / Re)

    def test_convergence_is_relative_change(self):
        geom = _fluid_geom()
        s    = geom.shape
        old  = 1.0 / Re
        new  = 2.0 / Re
        conv, _ = compute_convergence(_vel(2.0, s), geom, old)
        assert conv == pytest.approx(abs(new - old) / new)

    def test_from_zero_gives_convergence_one(self):
        # |new - 0| / new = 1
        geom = _fluid_geom()
        s    = geom.shape
        conv, _ = compute_convergence(_vel(2.0, s), geom, 0.0)
        assert conv == pytest.approx(1.0)

    def test_velocity_doubled_doubles_permeability(self):
        geom = _fluid_geom()
        s    = geom.shape
        _, p1 = compute_convergence(_vel(1.0, s), geom, 0.0)
        _, p2 = compute_convergence(_vel(2.0, s), geom, 0.0)
        assert p2 == pytest.approx(2.0 * p1)

    def test_all_solid_returns_unchanged_permeability(self):
        geom = make_proc_geom(np.ones((6, 6, 6), dtype=bool))
        s    = geom.shape
        old  = 0.42
        conv, new_perm = compute_convergence(_vel(1.0, s), geom, old)
        assert conv     == pytest.approx(0.0)
        assert new_perm == pytest.approx(old)

    def test_convergence_approaches_zero_over_iterations(self):
        # Simulate successive iterations with slowly changing velocity.
        geom = _fluid_geom()
        s    = geom.shape
        perm = 0.0
        vel  = 1.0
        for _ in range(10):
            conv, perm = compute_convergence(_vel(vel, s), geom, perm)
            vel *= 1.001   # tiny change each "iteration"
        # After many iterations with tiny changes convergence should be small.
        assert conv < 0.01
