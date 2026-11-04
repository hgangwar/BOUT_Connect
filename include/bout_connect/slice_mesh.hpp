#pragma once

#include "bout_connect/bout_interface.hpp"

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <random>

namespace bout_connect
{

struct SliceMeshOptions
{
  double jitter = 0.0; // > 0: each interior vertex moves by U(-jitter, jitter)
  double shift = 0.0;  // d: every interior vertex moves by (+d, +d)
  unsigned seed = 20260101;
};

inline Omega_h::Mesh build_slice_mesh(const Omega_h::CommPtr& world,
                                      const SliceMeshOptions& o = {},
                                      int* nmoved = nullptr)
{
  Omega_h::Mesh mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, double(NXI - 1), double(NZ), 0.0,
                       2 * (NXI - 1), 2 * NZ, 0, false);
  const int nv = mesh.nverts();
  Omega_h::HostRead<Omega_h::Real> hc(mesh.coords());
  Omega_h::HostWrite<Omega_h::Real> nc(2 * nv);
  std::mt19937 gen(o.seed);
  std::uniform_real_distribution<double> u(o.jitter > 0 ? -o.jitter : -1.0,
                                           o.jitter > 0 ? o.jitter : 1.0);
  int moved = 0;
  for (int v = 0; v < nv; ++v) {
    const double x = hc[2 * v], y = hc[2 * v + 1];
    const bool bnd =
      x < 1e-9 || x > (NXI - 1) - 1e-9 || y < 1e-9 || y > NZ - 1e-9;
    const double dx = bnd ? 0.0 : ((o.jitter == 0.0 ? 0.0 : u(gen)) + o.shift),
                 dy = bnd ? 0.0 : ((o.jitter == 0.0 ? 0.0 : u(gen)) + o.shift);
    moved += !bnd;
    nc[2 * v] = IX0 + x + dx; // index space: ix = IX0 + x
    nc[2 * v + 1] = y + dy;   // iz
  }
  mesh.set_coords(Omega_h::Reals(nc.write()));
  if (nmoved)
    *nmoved = moved;
  return mesh;
}

} // namespace bout_connect
