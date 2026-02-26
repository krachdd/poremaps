# Stokes Solver Benchmarks

This directory contains the three benchmark cases used to validate POREMAPS against known analytical solutions.

## Folder structure

```
simulations/
└── 99_benchmarks/
    ├── 01_poiseuille/      # Hagen-Poiseuille flow in a periodic tube
    ├── 02_channel_flow/    # Pressure-driven flow in a rectangular channel
    └── 03_spherepackings/  # Isotropic sphere-packing arrays (BCC, FCC)
```

Each subfolder contains the geometry (`.raw`), input file (`.inp`), and — for the Poiseuille case — pre-computed output fields and a log file.

## How to run

Build the binary first (see [../README.md](../README.md)), then run a benchmark from its subfolder, e.g.:

```shell
cd simulations/99_benchmarks/01_poiseuille
mpirun -np 4 ../../../bin/POREMAPS input_ptube_54_54_50_vs_2e-05.inp
```

---

## Benchmark 1 — Hagen-Poiseuille Flow in a Tube (`01_poiseuille/`)

**Geometry:** circular tube embedded in a 54 × 54 × 50 voxel periodic domain (voxel size $\Delta x = 2 \times 10^{-5}$ m).

**Analytical reference:** for a circular tube of radius $R$ under an axial pressure gradient $\Delta p / L$ and dynamic viscosity $\mu$, the Hagen-Poiseuille velocity profile is:
$$v_3(r) = \frac{\Delta p}{4L\mu} (R^2 - r^2)$$
See Batchelor, 1967.

**What to check:** compare the computed $v_3$ field (from `velz_*.raw`) to the parabolic profile above, evaluated at the radial distance $r$ of each fluid voxel from the tube centre.

**Included output:** pre-computed velocity, pressure, voxel-neighbourhood, and domain-decomposition fields are included alongside the log file.

---

## Benchmark 2 — Rectangular Channel Flow (`02_channel_flow/`)

**Geometry:** full rectangular cross-section channel, 24 × 14 × 50 voxels (voxel size $\Delta x = 5 \times 10^{-5}$ m; physical cross-section 1.2 mm × 0.7 mm).

**Analytical reference:** the volumetric flow rate $Q$ through a rectangular duct has a closed-form Fourier-series solution that depends on the aspect ratio of the cross-section and the applied pressure gradient. See White and Majdalani, 2006.

**What to check:** integrate the computed $v_3$ field (`velz_*.raw`) over the cross-section and compare the resulting flow rate $Q$ to the analytical value.

---

## Benchmark 3 — Sphere Packings (`03_spherepackings/`)

**Geometries:** two ordered sphere-packing arrays at porosity $\phi = 0.3$ (solid volume fraction 0.7), each 100 × 100 × 100 voxels (voxel size $\Delta x = 1 \times 10^{-5}$ m):

| Input file | Lattice | Sphere diameter $D$ |
|---|---|---|
| `input_SPHERE3D_bcc_…` | Body-centred cubic (BCC) | 8.745 × 10⁻⁴ m |
| `input_SPHERE3D_fcc_…` | Face-centred cubic (FCC) | 6.940 × 10⁻⁴ m |

**Analytical reference:** the Kozeny-Carman equation relates the intrinsic permeability $k$ to sphere diameter and porosity:
$$k^{KC} = \frac{D^2}{c_{KC}} \frac{\phi^3}{(1-\phi)^2}$$
where $c_{KC}$ is the Kozeny-Carman constant. See Carman, 1997 and Kozeny, 1927.

**What to check:** compare the whole-domain permeability $k_{33}$ (column `wk33` in the log file) to $k^{KC}$ computed from the geometry parameters above.
