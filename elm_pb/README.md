# elm_pb test case (verbatim copy)

Copied unchanged from BOUT-dev `examples/elm-pb` (BOUT++ is LGPL-3.0):

- `elm_pb.cxx` - the physics model, compiled into `bc_elmpb`; not clang-formatted on purpose.
- `cbm18_dens8.grid_nx68ny64.nc` - the CBM18 tokamak grid (the example symlinks it from `tokamak-grids`).
- `../data/BOUT.inp` - the input used by the tests (grid path relative to `build/` and `data/`; run from `build/` with `-d ../data`)
