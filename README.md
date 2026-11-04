# BOUT_Connect

Tests that bring a live BOUT++ field into PCMS **in memory**, in one process.

The BOUT++ `elm-pb` tokamak model is run unmodified. At the last output step the
solver is paused inside a monitor, PCMS reads `P` from BOUT++'s memory, and the
field is transferred and checked while BOUT++ waits.

Coordinates follow three conventions: index space `(ix, iy, iz)`, Cartesian
`(x, y, z)` and cylindrical `(R, phi, Z)`. The field is stored in index space
(`ix` 2..65, `iy` 0..63, `iz` 0..15, node 16 is the periodic copy of 0).

## Tests

| Target | What it does |
|---|---|
| `coexist` | Starts Kokkos and BOUT++ in one process and runs the elm-pb solver. A monitor does a small PCMS 3D transfer (3x3x3 to 5x5x5 nodes, `f = 1 + 2a + 3b + 4c`), which trilinear interpolation must reproduce to rounding. |
| `transfer` | Captures `P` at the pause, stores it as a PCMS 3D field in index space, transfers it to a 2x finer grid of the same box (532,257 nodes), and checks the result. |

`transfer` logs these steps (also written to `transfer_result.txt`):

- `[capture]` P copied out of BOUT++ memory.
- `[identity]` identity transfer against the BOUT++ buffer.
- `[2x]` 2x transfer: nodes that coincide with source nodes against BOUT++.
- Interpolator result against a PCMS point evaluation of the source at every target node.
- Value bounds of source and target.
- BOUT++'s own `XZBilinear` against PCMS at fractional-`ix`/`iz` points (20,000 target nodes, and 1,024 random points in the `iy = 32` row).
- `[norm]` nodal RMS and max error, against the PCMS source evaluation and against BOUT++.

There are no hand-written reference formulas: every comparison is against a
value BOUT++ computed or a value PCMS evaluated.

## Layout

| Path | Contents |
|---|---|
| `CMakeLists.txt` | Builds `bc_elmpb` (BOUT++ half), `coexist`, `transfer`. |
| `elm_pb/` | Verbatim copy of the BOUT-dev `examples/elm-pb` model and its grid (LGPL-3.0). See `elm_pb/README.md`. |
| `data/BOUT.inp` | The input used by the tests. The grid path is relative to `data/`. |
| `data/*.csv` | Output of the `transfer` run, tracked: `transfer_slice.csv`, `transfer_fine_slice.csv`, `transfer_targets.csv`, `transfer_random_targets.csv`. |
| `include/bout_connect/bout_interface.hpp` | The interface between the two halves: grid constants, `BoutSnapshot` (P copied out of BOUT++), the `run_elmpb*` entry points. Includes neither BOUT++ nor PCMS. |
| `include/bout_connect/bout_runner.hpp` | BOUT++-only: pauses the solver in a monitor, `capture()`, and the wrapper around BOUT++'s `XZBilinear`. |
| `include/bout_connect/pcms_util.hpp` | PCMS helpers: node coordinates, setting and getting DOFs, point evaluation, index-space function spaces. |
| `include/bout_connect/index_space.hpp` | The index-space coordinate system. |
| `include/bout_connect/slice_mesh.hpp` | Omega_h triangular-mesh builder; not used by `coexist` or `transfer`. |
| `include/bout_connect/report.hpp` | Logging and CSV reading. |
| `src/elmpb.cpp` | The BOUT++ half: includes the elm-pb source and implements the entry points. |
| `tests/` | `coexist.cpp` and `transfer.cpp`. |
| `plot_scripts/` | `plot_transfer.py` draws the verification figure from `data/*.csv`; `plots_dir.py` is its helper. Figures go to `plot_scripts/plots/`. |

### Why two halves

BOUT++ and PCMS cannot be compiled in the same translation unit: with a CUDA
Kokkos build, anything linking `pcms::core` is compiled by `nvcc`, and BOUT++'s
bundled `fmt` does not compile under `nvcc`. So `src/elmpb.cpp` (host compiler,
static library `bc_elmpb`) holds all BOUT++ code, each test (`nvcc`) holds all
PCMS code, and they talk only through plain C++ types and `std::function`
callbacks in `bout_interface.hpp`.

## Requirements

- A built BOUT++ (provides the `bout++` CMake package), 5.x.
- An installed PCMS that provides the `pcms` CMake package (Kokkos-enabled; the tests were run with the PETSc and CUDA build).
- MPI compilers, and `nvcc` on `PATH` if PCMS was built with CUDA Kokkos.
- Python with `numpy` and `matplotlib` for the plots.

## Build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=<mpicxx> -DCMAKE_C_COMPILER=<mpicc> \
      -Dbout++_DIR=<bout build or install dir> \
      -Dpcms_DIR=<pcms install>/lib64/cmake/pcms
cmake --build build -j 8
```

`CMakeLists.txt` sets `CMAKE_CXX20_STANDARD_COMPILE_OPTION` to `-std=c++20`,
because `nvcc_wrapper` otherwise downgrades `-std=c++2a` to C++17 and a PCMS
header needs C++20. To skip a test, comment out its `add_executable` and
`target_link_libraries` lines in `CMakeLists.txt`.

## Run

BOUT++ reads its input and grid relative to the working directory, so run from
`data/` with `-d .`:

```
cd data
../build/coexist  -d .
../build/transfer -d . solver:nout=1
```

`solver:nout=1` pauses at output 1 (t = 1). The committed CSVs come from that
run; without it the input's default `nout = 20` pauses at t = 20, where `P` is
larger and the CSVs differ. The BOUT++ outputs (`BOUT.dmp`, `BOUT.restart`,
logs) are ignored by git.

## Plot

```
cd plot_scripts
python plot_transfer.py
```

It reads the CSVs from `<package>/data` (or `$BC_DATA`) and writes
`plot_scripts/plots/transfer_verification.png`.

## Expected result

With `solver:nout=1`: `transfer` exits 0; the identity and 2x transfers match
BOUT++ exactly (max diff 0); the Interpolator matches the PCMS point evaluation
at all 532,257 target nodes (max diff 0); BOUT++'s `XZBilinear` agrees with PCMS
to about 1e-16 of max|P|. `coexist` prints `max |target - f| = 0.000e+00`.
