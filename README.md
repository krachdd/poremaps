[![POREMAPS](doc/poremaps_logo_grey.png)](https://www.mib.uni-stuttgart.de/cont/)

**POREMAPS** is a Finite Difference based <ins>POR</ins>ous <ins>M</ins>edia <ins>A</ins>nisotropic <ins>P</ins>ermeability <ins>S</ins>olver for Stokes flow.

[![Identifier](https://img.shields.io/badge/doi-10.18419%2Fdarus--3676-d45815)](https://doi.org/10.18419/darus-3676)
[![Identifier](https://img.shields.io/badge/Publication-blue)](https://doi.org/10.69631/ipj.v2i1nr39)

Given a binarized 3-D pore geometry (`.raw` file), POREMAPS solves the Stokes equations on a Cartesian finite-difference grid using MPI domain decomposition and computes the anisotropic intrinsic permeability tensor of the sample.

---

## Contents

- [How to build](#how-to-build)
- [How to run](#how-to-run)
- [Testing](#testing)
- [Boundary conditions](#boundary-conditions)
- [Input data](#input-data)
- [Output files](#output-files)
- [Aspects of parallelization](#aspects-of-parallelization)
- [Compute the permeability tensor](#compute-the-permeability-tensor)
- [License](#license)
- [How to cite](#how-to-cite)

---

## How to build

POREMAPS requires an MPI C++ compiler (`mpiCC`). We recommend building a local version of [OpenMPI](https://www.open-mpi.org/). Tested versions:

| Platform | MPI implementation |
|---|---|
| Linux | `openmpi-3.1.5`, `openmpi-4.1.5`, `openmpi-5.0.1` |
| Windows 10 | `MS-MPI v10.1.1` |

**Linux** — if `mpiCC` is on your `PATH` (standard OpenMPI installation), just run:
```shell
cd PATH/TO/POREMAPS/src/
make
```

If your MPI compiler is installed at a non-standard location, pass it on the command line:
```shell
make MPICXX=/path/to/your/mpiCC
```

A debug build (no optimisation, debug symbols) is also available:
```shell
make debug
```

The executable `POREMAPS` is placed in the `bin/` directory.

**Windows** — create a Visual Studio (v16.8.3) project from the source files and link MS-MPI (v10.1.1) via the project properties.

---

## How to run

```shell
mpirun -np 8 ./bin/POREMAPS my_inputfile.inp
```

The input file must follow the exact format shown in `input_template.inp`. All parameters are read in order; the field names are ignored by the parser.

```text
dom_decomposition 0 0 0
boundary_method 0
geometry_file_name ptube_54_54_50_vs_2e-05.raw
size_x_y_z  54 54 50
voxel_size  2e-05
max_iter    100000001
it_eval    100
it_write   100
log_file_name permeability_ptube_54_54_50_vs_2e-05.log
solving_algorithm 2
eps 1e-06
porosity 1.0
dom_interest 0 0 0 0 0 0
write_output 1 1 0 0
```

### Parameter reference

| Parameter | Value(s) | Description |
|---|---|---|
| `dom_decomposition` | `int int int` | Number of MPI ranks in $\mathbf{e}_1$, $\mathbf{e}_2$, $\mathbf{e}_3$. Set to `0 0 0` to let MPI choose automatically. |
| `boundary_method` | `0`–`4` | Physical boundary setup; see [Boundary conditions](#boundary-conditions). |
| `geometry_file_name` | `string` | Path to the input geometry file (`.raw`). |
| `size_x_y_z` | `int int int` | Global domain size in voxels $(n_x,\ n_y,\ n_z)$. |
| `voxel_size` | `float` | Voxel edge length in metres. |
| `max_iter` | `int` | Maximum number of pseudo-time iterations. |
| `it_eval` | `int` | Evaluate the convergence criterion every `it_eval` iterations. |
| `it_write` | `int` | Write a log entry every `it_write` iterations (must be a multiple of `it_eval`). |
| `log_file_name` | `string` | Path for the output log file. |
| `solving_algorithm` | `1`, `2`, or `3` | `1` = Jacobi, `2` = Gauss–Seidel, `3` = SOR. |
| `eps` | `float` | Convergence threshold; iteration stops when the relative permeability change falls below this value. |
| `porosity` | `float` | Porosity of the sample. Set to `-1.0` to have the solver compute it from the geometry (not the effective porosity). Specify explicitly when using solid frames or sub-domain evaluation. |
| `dom_interest` | `x_min x_max y_min y_max z_min z_max` | Voxel index bounds of a sub-domain of interest. The solver reports additional permeability values for this region. Set all entries to `0` to skip sub-domain evaluation. |
| `write_output` | `int int int int` | Controls which field files are written: velocity (`velx/y/z`), pressure, voxel neighbourhood, domain decomposition — in that order. `1` = write, `0` = skip. |

---

## Testing

A Python test suite covering the core numerical kernels lives in the `tests/` directory. The tests verify the finite-difference stencils, domain decomposition arithmetic, boundary condition assignment, and permeability and convergence formulae via reference implementations that mirror the C++ source without requiring MPI.

**Dependencies** — `numpy` and `pytest`:
```shell
pip install numpy pytest
```

**Run all unit tests** from the repository root:
```shell
pytest tests/
```

**Run with verbose output:**
```shell
pytest tests/ -v
```

**Integration tests** invoke the compiled `POREMAPS` binary on small generated geometries. They require `mpirun` on `PATH` and the binary to be built. They are opt-in:
```shell
pytest tests/ --run-integration
```

| Module | Functions tested | Source file |
|---|---|---|
| `test_solver.py` | `split_number`, `get_2nd_derivation` | `solver.cc` |
| `test_geometry.py` | `eval_geometry`, `get_dom_limits`, `get_proc_porosity` | `geometry.cc` |
| `test_evaluation.py` | `compute_permeability`, `compute_convergence` | `evaluation.cc` |
| `test_initialize.py` | `determine_comm_pattern`, `set_halos_initial` | `initialize.cc` |
| `test_integration.py` | Full binary run via `mpirun` | — (opt-in) |

---

## Boundary conditions

The main pressure gradient and flow direction is always $\mathbf{e}_3$ (z). Choose `boundary_method` according to the physical setup:

| Value | z-direction | x- and y-directions |
|:-----:|---|---|
| `0` | periodic | periodic |
| `1` | periodic | slip walls |
| `2` | periodic | no-slip walls |
| `3` | non-periodic (pressure-driven) | slip walls |
| `4` | non-periodic (pressure-driven) | no-slip walls |

---

## Input data

POREMAPS reads 8-bit binary (`.raw`) geometry files. Each byte encodes one voxel: `0` = fluid, `1` = solid. There is no file header — dimensions must be provided via `size_x_y_z` in the input file.

**Byte order:** data must be stored with x varying fastest (Fortran index order). When preparing input files with NumPy, use a `(nz, ny, nx)` shaped array written in C order, which gives x-fastest layout:

```python
import numpy as np
geom = np.zeros((nz, ny, nx), dtype=np.uint8)  # 0 = fluid, 1 = solid
geom.tofile("geometry.raw")  # default C order → x varies fastest
```

**Connectivity:** the geometry must contain at least one connected flow path from the inlet to the outlet face (face-to-face connectivity). This can be verified with:

```python
from scipy.ndimage import label
labeled, n = label(geom == 0)  # label connected fluid regions
```

Preprocessing (mirroring, padding, solid frame addition) must be done before passing the file to POREMAPS. We recommend NumPy for this.

To visualise `.raw` output fields, use the included `fields2vtu.py` script to convert them to `.vtu` format for [ParaView](https://www.paraview.org/).

---

## Output files

### Field files (`.raw`)

Written when the corresponding `write_output` flag is set to `1`:

| File pattern | Data type | Content |
|---|---|---|
| `velx_*.raw`, `vely_*.raw`, `velz_*.raw` | `double` | Final velocity field ($\mathbf{e}_1$, $\mathbf{e}_2$, $\mathbf{e}_3$ components) |
| `press_*.raw` | `double` | Final pressure field |
| `voxel_neighborhood_*.raw` | `int` | Voxel neighbourhood case codes |
| `domain_decomp_*.raw` | `int` | MPI domain decomposition map |

Use `fields2vtu.py` to convert any of these to `.vtu` for visualisation in [ParaView](https://www.paraview.org/).

### Log file (`.log`)

One row per write event; columns in order:

| Column | Unit | Description |
|---|---|---|
| `iteration` | — | Iteration number |
| `conv` | — | Relative permeability change (convergence criterion) |
| `TPS` | 1/s | Time steps computed per wall-clock second |
| `wmax_velz` | m/s | Maximum fluid-voxel z-velocity (whole domain) |
| `wmean_velz` | m/s | Mean fluid-voxel z-velocity (whole domain) |
| `k13`, `k23`, `k33` | m² | Permeability components for the `dom_interest` sub-domain (`0.0` if not set) |
| `wk13`, `wk23`, `wk33` | m² | Permeability components for the whole domain |

---

## Aspects of parallelization

A useful rule of thumb is to target **50³–100³ voxels per MPI rank** (smaller sub-domains for high-porosity samples, larger for low-porosity ones). Larger sub-domains are technically possible and may be preferable when memory bandwidth is the bottleneck.

To assess parallel efficiency for a given problem, run a short simulation (a few hundred iterations) and inspect the TPS value in the log. Experimenting with manual `dom_decomposition` settings versus the automatic (`0 0 0`) decomposition can also be worthwhile for non-cubic domains. Core counts used for larger simulations (> 1500³ voxels) are reported in the paper below.

---

## Compute the permeability tensor

The pressure gradient is always applied along $\mathbf{e}_3$ (z). To obtain off-diagonal tensor entries, rotate the domain and run the solver once per orientation. Each run fills one column of the permeability tensor:

| `numpy.transpose` call | Output column |
|---|---|
| *(no transpose — original)* | $k_{13},\ k_{23},\ k_{33}$ |
| `np.transpose(domain, (2, 0, 1))` | $k_{32},\ k_{12},\ k_{22}$ |
| `np.transpose(domain, (1, 2, 0))` | $k_{21},\ k_{31},\ k_{11}$ |

Remember to rotate the result vectors by the same permutation to map log-file entries back to physical tensor indices.

---

## License

POREMAPS is released under the [MIT License](https://opensource.org/license/mit/). See `LICENSE.md` in the repository root for the full text.

---

## How to cite

Please cite **specific releases** via [DaRUS](https://doi.org/10.18419/darus-3676).

If you use POREMAPS in a scientific publication, please also cite:

```bib
@article{Krach2025a,
    author  = {Krach, David and Ruf, Matthias and Steeb, Holger},
    title   = {{POREMAPS}: A finite difference based Porous Media Anisotropic Permeability Solver for Stokes flow},
    DOI     = {10.69631/ipj.v2i1nr39},
    journal = {InterPore Journal},
    pages   = {IPJ260225--7},
    volume  = {2},
    number  = {1},
    month   = {Feb.},
    year    = {2025},
    place   = {De Bilt, The Netherlands}
}
```

```bib
@data{Krach2024a,
    author    = {Krach, David and Ruf, Matthias and Steeb, Holger},
    publisher = {DaRUS},
    title     = {{POREMAPS 1.0.0: Code, Benchmarks, Applications}},
    year      = {2024},
    version   = {V1},
    doi       = {10.18419/darus-3676},
    url       = {https://doi.org/10.18419/darus-3676}
}
```

---

## Solver is used in the following publications

[![Identifier](https://img.shields.io/badge/Publication_ADWR_Krach_et.al._(2025)-blue)](https://doi.org/10.1016/j.advwatres.2024.104860)

```bib
@article{Krach2025b,
    title   = {A novel geometry-informed drag term formulation for pseudo-3D Stokes simulations with varying apertures},
    journal = {Advances in Water Resources},
    volume  = {195},
    year    = {2025},
    doi     = {https://doi.org/10.1016/j.advwatres.2024.104860},
    author  = {David Krach and Felix Weinhardt and Mingfeng Wang and Martin Schneider and Holger Class and Holger Steeb},
    keywords = {Porous media, Stokes flow, Biomineralization, Microfluidics, Image-based simulations, Computational efficiency versus accuracy}
}
```

---

## Developer

- [David Krach](https://www.mib.uni-stuttgart.de/institute/team/Krach/) — [david.krach@mib.uni-stuttgart.de](mailto:david.krach@mib.uni-stuttgart.de)
- [Matthias Ruf](https://www.mib.uni-stuttgart.de/institute/team/Ruf-00001/) — [matthias.ruf@mib.uni-stuttgart.de](mailto:matthias.ruf@mib.uni-stuttgart.de)

## Contact

- [Software Support — Institute of Applied Mechanics](mailto:software@mib.uni-stuttgart.de)
- [Data Support — Institute of Applied Mechanics](mailto:data@mib.uni-stuttgart.de)
