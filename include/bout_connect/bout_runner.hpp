// BOUT++-only (host compiler): runs a physics model, pauses it in a monitor,
// and copies P out into a BoutSnapshot.
#pragma once

#include "bout_connect/bout_interface.hpp"
#include "bout_connect/report.hpp"

#include <bout/bout.hxx>
#include <bout/interpolation_xz.hxx>
#include <bout/monitor.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/solver.hxx>

#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace bout_connect
{
class PauseMonitor : public Monitor
{
public:
  using Work = std::function<void(Solver*, BoutReal, int, int)>;
  PauseMonitor(Work work, bool& ran) : work_(std::move(work)), ran_(ran) {}
  int call(Solver* solver, BoutReal t, int iter, int nout) override
  {
    if (iter == 0 || iter != nout)
      return 0;
    ran_ = true;
    try {
      work_(solver, t, iter, nout);
    } catch (const std::exception& e) {
      say("EXCEPTION inside the pause: %s\n", e.what());
    }
    return 0;
  }

private:
  Work work_;
  bool& ran_;
};

// Runs the BOUT++ model `Model` with the monitor
template <class Model>
int run_model(Monitor& mon)
{
  try {
    auto model = bout::utils::make_unique<Model>();
    auto solver = Solver::create();
    solver->setModel(model.get());
    auto bout_monitor = bout::utils::make_unique<BoutMonitor>();
    solver->addMonitor(bout_monitor.get(), Solver::BACK);
    solver->addMonitor(&mon, Solver::BACK);
    solver->solve();
  } catch (const std::exception& e) {
    say("EXCEPTION in the BOUT++ part: %s\n", e.what());
    return 2;
  }
  return 0;
}

// run_model with a PauseMonitor
template <class Model>
int run_paused(PauseMonitor::Work work, bool& ran)
{
  PauseMonitor mon(std::move(work), ran);
  return run_model<Model>(mon);
}

// P (and the geometry tables) copied out of BOUT++ memory.
struct Captured
{
  BoutSnapshot data;
  Field3D field;
};

inline Captured capture(Solver* solver)
{
  Mesh* mesh = bout::globals::mesh;
  Options o;
  solver->outputVars(o, false);
  Captured c;
  c.field = o["P"].as<Field3D>(Field3D{mesh});
  BoutSnapshot& d = c.data;
  for (const char* key : {"direction_y", "direction_z", "cell_location"})
    if (o["P"].hasAttribute(key)) {
      const std::string v = o["P"].attributes.at(key);
      d.attrs += std::string(key) + "=" + v + " ";
    }
  d.pb.resize(size_t(NXI) * NY * NZ);
  for (int ix = IX0; ix < IX0 + NXI; ++ix)
    for (int iy = 0; iy < NY; ++iy)
      for (int iz = 0; iz < NZ; ++iz) {
        const double v = c.field(ix, iy + mesh->ystart, iz);
        d.pb[size_t(iz) + NZ * (size_t(iy) + size_t(NY) * (ix - IX0))] = v;
        d.pmax = std::max(d.pmax, std::fabs(v));
        d.pmin = std::min(d.pmin, v);
        d.pmaxv = std::max(d.pmaxv, v);
      }
  Field2D Rxy{mesh}, Zxy{mesh};
  mesh->get(Rxy, "Rxy");
  mesh->get(Zxy, "Zxy");
  d.dz = getUniform(mesh->getCoordinates()->dz());
  d.R.resize(size_t(NXI) * NY);
  d.Z.resize(size_t(NXI) * NY);
  for (int ix = IX0; ix < IX0 + NXI; ++ix)
    for (int iy = 0; iy < NY; ++iy) {
      d.R[size_t(ix - IX0) + size_t(NXI) * iy] = Rxy(ix, iy + mesh->ystart);
      d.Z[size_t(ix - IX0) + size_t(NXI) * iy] = Zxy(ix, iy + mesh->ystart);
    }
  return c;
}

// BOUT++'s XZBilinear of P at points (ix, iy_file, iz)
inline std::vector<double> xz_bilinear(
  const Field3D& P, const std::vector<std::array<double, 3>>& q)
{
  Mesh* mesh = bout::globals::mesh;
  Field3D dx{0.0}, dz{0.0};
  for (int jx = 0; jx < mesh->LocalNx; ++jx)
    for (int jy = 0; jy < mesh->LocalNy; ++jy)
      for (int jz = 0; jz < mesh->LocalNz; ++jz) {
        dx(jx, jy, jz) = jx;
        dz(jx, jy, jz) = jz;
      }
  struct Hold
  {
    int ix, iy, iz;
  };
  std::vector<int> rowcnt(NY, 0);
  std::vector<Hold> holds(q.size());
  for (size_t k = 0; k < q.size(); ++k) {
    const int iy = int(std::lround(q[k][1]));
    const int cnt = rowcnt[iy]++;
    if (cnt >= NXI * NZ)
      throw std::runtime_error(
        "xz_bilinear: more than NXI * NZ points in one y row");
    holds[k] = {IX0 + cnt % NXI, iy + mesh->ystart, (cnt / NXI) % NZ};
    dx(holds[k].ix, holds[k].iy, holds[k].iz) =
      q[k][0]; // absolute fractional ix
    dz(holds[k].ix, holds[k].iy, holds[k].iz) = q[k][2]; // and iz of the point
  }
  XZBilinear bil(0, mesh);
  bil.calcWeights(dx, dz);
  const Field3D fb = bil.interpolate(P);
  std::vector<double> out(q.size());
  for (size_t k = 0; k < q.size(); ++k)
    out[k] = fb(holds[k].ix, holds[k].iy, holds[k].iz);
  return out;
}
} // namespace bout_connect
