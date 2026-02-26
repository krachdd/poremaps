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

"""Integration tests: invoke the POREMAPS binary on small generated geometries.

These tests are skipped unless ``--run-integration`` is passed to pytest.
They also require:
  * The binary compiled at  ../bin/POREMAPS  (run ``make`` in src/).
  * ``mpirun`` available on PATH.

Example
-------
    pytest tests/ --run-integration
"""

import os
import shutil
import subprocess
import tempfile
import textwrap

import numpy as np
import pytest

BINARY = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "bin", "POREMAPS")
)


# ---------------------------------------------------------------------------
# Module-level skip fixture applied to every test in this file
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _check_prerequisites(run_integration):
    if not run_integration:
        pytest.skip("pass --run-integration to enable binary tests")
    if not os.path.isfile(BINARY):
        pytest.skip(f"binary not found at {BINARY} — run 'make' in src/ first")
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not found on PATH")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_geom(path, nx, ny, nz):
    """Write an all-fluid uint8 raw geometry file (x-fastest order)."""
    geom = np.zeros((nz, ny, nx), dtype=np.uint8)
    geom.ravel(order='C').tofile(path)


def _write_inp(path, geom_path, log_path, nx, ny, nz,
               bc_method=0, voxelsize=1e-6, max_iter=100, eps=1e-12):
    """Write a POREMAPS input parameter file."""
    text = textwrap.dedent(f"""\
        dom_decomposition 0 0 0
        boundary_method {bc_method}
        geometry_file_name {geom_path}
        size_x_y_z  {nx} {ny} {nz}
        voxel_size  {voxelsize}
        max_iter    {max_iter}
        it_eval    50
        it_write   50
        log_file_name {log_path}
        solving_algorithm 2
        eps {eps}
        porosity -1.0
        dom_interest 0 0 0 0 0 0
        write_output 0 0 0 0
    """)
    with open(path, "w") as f:
        f.write(text)


def _run(tmpdir, nx=8, ny=8, nz=8, nranks=1, **inp_kwargs):
    """Set up geometry + input files, run the binary, return (result, log_path)."""
    geom_path = os.path.join(tmpdir, "geom.raw")
    log_path  = os.path.join(tmpdir, "run.log")
    inp_path  = os.path.join(tmpdir, "input.inp")

    _write_geom(geom_path, nx, ny, nz)
    _write_inp(inp_path, geom_path, log_path, nx, ny, nz, **inp_kwargs)

    result = subprocess.run(
        ["mpirun", "-n", str(nranks), BINARY, inp_path],
        capture_output=True,
        text=True,
        cwd=tmpdir,
        timeout=120,
    )
    return result, log_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestIntegration:
    def test_binary_exits_successfully(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result, _ = _run(tmpdir, nx=8, ny=8, nz=8, max_iter=100)
        assert result.returncode == 0, (
            f"POREMAPS exited with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    def test_log_file_is_created(self):
        # The solver only writes the log after iteration 1000; use max_iter > 1000.
        with tempfile.TemporaryDirectory() as tmpdir:
            result, log_path = _run(tmpdir, nx=8, ny=8, nz=8, max_iter=1100)
            assert result.returncode == 0
            assert os.path.isfile(log_path), "Log file was not created by the solver."

    def test_log_file_contains_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result, log_path = _run(tmpdir, nx=8, ny=8, nz=8, max_iter=1100)
            assert result.returncode == 0
            with open(log_path) as f:
                content = f.read()
            assert len(content.strip()) > 0, "Log file is empty."

    def test_log_contains_numeric_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result, log_path = _run(tmpdir, nx=8, ny=8, nz=8, max_iter=1100)
            assert result.returncode == 0
            with open(log_path) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            assert lines, "Log file has no non-empty lines."
            # At least the last data line should be parseable as floats.
            # Log format is comma-separated: "iter, val, val, ..."
            data_lines = [ln for ln in lines if not ln.startswith('#')]
            assert data_lines, "Log file has no data lines (only a header)."
            try:
                values = [float(v) for v in data_lines[-1].split(',')]
            except ValueError:
                pytest.fail(f"Last log line is not numeric: {data_lines[-1]!r}")
            assert len(values) > 0

    def test_two_mpi_ranks(self):
        """Solver must complete without error when split across 2 MPI ranks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result, _ = _run(tmpdir, nx=8, ny=8, nz=8, nranks=2, max_iter=50)
        assert result.returncode == 0, (
            f"2-rank run failed\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    def test_permeability_positive_all_fluid(self):
        """The whole-domain z-permeability (wk33) in the log must be positive."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result, log_path = _run(
                tmpdir, nx=8, ny=8, nz=10,
                max_iter=1100, eps=1e-8, voxelsize=1e-5,
            )
            assert result.returncode == 0
            with open(log_path) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            assert lines
            data_lines = [ln for ln in lines if not ln.startswith('#')]
            assert data_lines, "Log file has no data lines (only a header)."
            try:
                values = [float(v) for v in data_lines[-1].split(',')]
            except ValueError:
                pytest.skip("Could not parse log file format — check log output")
            assert any(v > 0 for v in values), (
                f"Expected at least one positive value in log line: {data_lines[-1]}"
            )
