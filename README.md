# BOUT_Connect

Brings the live BOUT++ `elm-pb` field into PCMS in memory, in one process. The
solver is paused in a monitor at the last output step while PCMS works.

## Tests

- `coexist`: Kokkos and BOUT++ in one process; a small PCMS 3D transfer inside the monitor.
- `transfer`: captures `P`, stores it as a PCMS 3D field in index space `(ix, iy, iz)`, transfers it to a 2x finer grid, and compares with BOUT++ values and PCMS evaluations.

## Layout

- `elm_pb/`: Exact BOUT-dev elm-pb model and grid.
- `data/`: `BOUT.inp` and the `transfer` CSV output.
- `include/bout_connect/`: `bout_interface.hpp` (BOUT++-free interface), `bout_runner.hpp` (BOUT++ side), `pcms_util.hpp`, `index_space.hpp`, `report.hpp`, `slice_mesh.hpp`.
- `src/elmpb.cpp`: BOUT++ half, built with the host compiler.
- `tests/`: `coexist.cpp`, `transfer.cpp` (built with nvcc, PCMS half).
- `plot_scripts/`: `plot_transfer.py`.

BOUT++ and PCMS live in separate translation units because BOUT++'s `fmt` does not compile under nvcc.

## Build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dbout++_DIR=<bout build dir> -Dpcms_DIR=<pcms install>/lib64/cmake/pcms
cmake --build build -j 8
```

## Run (from `build/`)

```
cd build
./coexist -d ../data
./transfer -d ../data nout=1
```

`-d ../data` points BOUT++ at `data/BOUT.inp`; the CSVs and `transfer_result.txt` are written to `data/`.

## Plot

```
python plot_scripts/plot_transfer.py   # reads data/*.csv
```
