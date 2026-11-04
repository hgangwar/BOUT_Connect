// coexist: BOUT++ and PCMS in ONE process.
#include "bout_connect/bout_interface.hpp"
#include "bout_connect/pcms_util.hpp"

#include <Kokkos_Core.hpp>
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/transfer/interpolator.h"
#include "pcms/utility/uniform_grid.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace bout_connect;

namespace
{
double f_lin(double a, double b, double c)
{ // the linear test function f = 1 + 2a + 3b + 4c
  return 1.0 + 2.0 * a + 3.0 * b + 4.0 * c;
}

struct Result
{
  bool done = false;
  double maxerr = 0;
  size_t ntgt = 0;
};
Result g_res;

Result pcms_trivial_transfer()
{
  Result r;
  pcms::UniformGrid<3> gs, gt;
  gs.bot_left = {0.0, 0.0, 0.0};
  gs.edge_length = {2.0, 2.0, 2.0};
  gs.divisions = {2, 2, 2}; // 3x3x3 nodes
  gt.bot_left = {0.0, 0.0, 0.0};
  gt.edge_length = {2.0, 2.0, 2.0};
  gt.divisions = {4, 4, 4}; // 5x5x5 nodes, same box
  const auto cs = pcms::csys::Cartesian::Create(3);
  auto src_space = pcms::LagrangeFunctionSpace::FromUniformGrid(gs, 1, cs, 1);
  auto tgt_space = pcms::LagrangeFunctionSpace::FromUniformGrid(gt, 1, cs, 1);
  auto src = src_space->CreateFunction<pcms::Real>();
  auto tgt = tgt_space->CreateFunction<pcms::Real>();

  // fill the source from the source space's own node coordinates
  const Pts<3> sn = node_coords<3>(*src_space);
  const size_t ns = sn.size();
  std::vector<double> sv(ns);
  for (size_t k = 0; k < ns; ++k)
    sv[k] = f_lin(sn[k][0], sn[k][1], sn[k][2]);
  set_dofs(src, sv);

  pcms::Interpolator<pcms::Real> interp(*src_space, *tgt_space);
  interp.Apply(src, tgt);

  const Pts<3> tn = node_coords<3>(*tgt_space);
  r.ntgt = tn.size();
  const auto tv = get_dofs(tgt, r.ntgt);
  for (size_t k = 0; k < r.ntgt; ++k)
    r.maxerr = std::max(r.maxerr,
                        std::fabs(tv[k] - f_lin(tn[k][0], tn[k][1], tn[k][2])));
  r.done = true;
  return r;
}

// Called by the BOUT++ solver
void on_monitor_call(double t, int iter, int nout)
{
  std::printf("\n[coexist monitor] called by the BOUT++ solver: t = %g, "
              "iteration %d of %d\n",
              t, iter, nout);
  if (g_res.done)
    return;
  try {
    g_res = pcms_trivial_transfer();
  } catch (const std::exception& e) {
    std::printf("[coexist monitor] PCMS EXCEPTION: %s\n", e.what());
    g_res.done = true;
    return;
  }
  std::printf("[coexist monitor] PCMS 3D transfer 3x3x3 -> 5x5x5 (%zu "
              "target nodes): max |target - f| = %.3e\n",
              g_res.ntgt, g_res.maxerr);
}
} // namespace

int main(int argc, char** argv)
{
  Kokkos::ScopeGuard kokkos(argc, argv);
  std::printf("Kokkos default execution space: %s;  concurrency %d\n",
              Kokkos::DefaultExecutionSpace::name(),
              int(Kokkos::DefaultExecutionSpace().concurrency()));
#ifdef _OPENMP
  std::printf(
    "OpenMP enabled in this translation unit: omp_get_max_threads() = %d\n",
    omp_get_max_threads());
#else
  std::printf("OpenMP not enabled in this translation unit\n");
#endif
  const int status = run_elmpb_monitor(argc, argv, on_monitor_call);
  std::printf("\nco-existence run finished: the monitor %s, BOUT++ exit status "
              "%d\n",
              g_res.done ? "ran the PCMS transfer" : "did NOT run", status);
  return g_res.done ? status : 1;
}
