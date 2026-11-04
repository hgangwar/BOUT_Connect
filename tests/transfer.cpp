// transfer: BOUT++ to PCMS in ONE process.  The real elm-pb run is paused
// (inside a monitor, between solver steps) at the last output step; at that
// pause PCMS reads P from BOUT++'s memory, stores it as a 3D field in INDEX
// space, transfers it to a 2x finer grid of the same box, and the result is
// checked against the source.

#include "bout_connect/bout_interface.hpp"
#include "bout_connect/pcms_util.hpp"
#include "bout_connect/report.hpp"

#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/value_view.hpp"
#include "pcms/transfer/interpolator.h"
#include "pcms/utility/uniform_grid.h"
#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

using namespace bout_connect;

namespace
{
void work(const BoutSnapshot& PB, double t, int iter, int nout)
{
  using clock = std::chrono::steady_clock;
  const auto t_pause = clock::now();
  say("\n==== BOUT++ paused at output %d of %d (t = %g); PCMS works now, "
      "inside the monitor ====\n",
      iter, nout, t);

  // ---------------- capture: P from BOUT++ memory
  const double pmax = PB.pmax, pmin = PB.pmin, pmaxv = PB.pmaxv;
  const std::string& attrs = PB.attrs;
  say(
    "[capture] P captured: interior %dx%dx%d, min %.6e max %.6e max|P| %.6e;  "
    "attributes: %s\n",
    NXI, NY, NZ, pmin, pmaxv, pmax, attrs.c_str());

  // ---------------- identity: PCMS source field + identity transfer
  auto t0 = clock::now();
  auto src_space = make_index_space(1);
  auto src = src_space->CreateFunction<pcms::Real>();
  const Pts<3> sn = node_coords<3>(*src_space);
  std::vector<double> sv(sn.size());
  for (size_t k = 0; k < sn.size(); ++k)
    sv[k] = PB(int(std::lround(sn[k][0])) - IX0, int(std::lround(sn[k][1])),
               int(std::lround(sn[k][2]))); // (ix, iy, iz)
  set_dofs(src, sv);
  auto same_space = make_index_space(1);
  auto same = same_space->CreateFunction<pcms::Real>();
  pcms::Interpolator<pcms::Real> ident(*src_space, *same_space);
  ident.Apply(src, same);
  const auto sb = get_dofs(same, sn.size());
  double d3 = 0;
  for (size_t k = 0; k < sn.size(); ++k)
    d3 = std::max(d3, std::fabs(sb[k] - sv[k]));
  say("[identity] source field: %zu nodes (ix 2..65, iy 0..63, iz 0..16);  "
      "identity transfer vs BOUT++ buffer: max diff = %.3e   [%.2f s]\n",
      sn.size(), d3, secs(t0));

  // ---------------- 2x: transfer
  t0 = clock::now();
  auto tgt_space = make_index_space(2);
  auto tgt = tgt_space->CreateFunction<pcms::Real>();
  pcms::Interpolator<pcms::Real> xfer(*src_space, *tgt_space);
  xfer.Apply(src, tgt);
  const Pts<3> tn = node_coords<3>(*tgt_space);
  const auto tv = get_dofs(tgt, tn.size());
  double d4 = 0;
  size_t ncoin = 0;
  for (size_t k = 0; k < tn.size(); ++k)
    if (is_int(tn[k][0]) && is_int(tn[k][1]) && is_int(tn[k][2])) {
      ++ncoin;
      d4 = std::max(d4, std::fabs(tv[k] - PB(int(std::lround(tn[k][0])) - IX0,
                                             int(std::lround(tn[k][1])),
                                             int(std::lround(tn[k][2])))));
    }
  say("[2x] 2x transfer: %zu target nodes (%d x %d x %d);  %zu coincide with "
      "source nodes, max diff to BOUT++ = %.3e   [%.2f s]\n",
      tn.size(), 2 * NZ + 1, 2 * (NY - 1) + 1, 2 * (NXI - 1) + 1, ncoin, d4,
      secs(t0));

  // PCMS point evaluation of the SOURCE field at target nodes
  t0 = clock::now();
  const auto ev_src =
    eval_points<3>(src_space, src, tn, IndexSpace3D::Create());
  double d5b = 0, tmin = 1e300, tmax = -1e300;
  for (size_t k = 0; k < tn.size(); ++k) {
    d5b = std::max(d5b, std::fabs(tv[k] - ev_src[k]));
    tmin = std::min(tmin, tv[k]);
    tmax = std::max(tmax, tv[k]);
  }
  say("Interpolator result vs PCMS point evaluation of the "
      "source at all %zu target nodes: max diff = %.3e  (relative to max|P|: "
      "%.3e)\n",
      tn.size(), d5b, d5b / pmax);
  say("bounds: target [%.6e, %.6e], source [%.6e, %.6e]\n", tmin, tmax, pmin,
      pmaxv);

  // BOUT++'s own XZBilinear at target nodes
  std::vector<size_t> cand;
  for (size_t k = 0; k < tn.size(); ++k)
    if (is_int(tn[k][1]) && tn[k][2] <= NZ - 1 + 1e-12 &&
        !(is_int(tn[k][0]) && is_int(tn[k][2])))
      cand.push_back(k); // iy integer, iz <= 15, ix or iz fractional
  std::mt19937 gen(777);
  std::shuffle(cand.begin(), cand.end(), gen);
  std::vector<int> rowcnt(NY, 0);
  std::vector<size_t> held;
  Pts<3> hq;
  for (size_t k : cand) {
    if (held.size() >= 20000)
      break;
    const int iy = int(std::lround(tn[k][1]));
    if (rowcnt[iy] >= NXI * NZ)
      continue;
    ++rowcnt[iy];
    held.push_back(k);
    hq.push_back(tn[k]);
  }
  const auto fb = PB.xz_bilinear(hq);
  double d5d = 0;
  for (size_t i = 0; i < held.size(); ++i)
    d5d = std::max(d5d, std::fabs(tv[held[i]] - fb[i]));
  say("BOUT++ XZBilinear at %zu target nodes (integer iy, "
      "fractional ix "
      "and/or iz): max diff = %.3e  (relative: %.3e)\n",
      held.size(), d5d, d5d / pmax);
  say("     (checks b, c, d: %.2f s)\n", secs(t0));

  // Data for the verification plot
  {
    const int ROW = 32;
    if (FILE* fp = std::fopen(data_path("transfer_slice.csv").c_str(),
                              "w")) { // the SOURCE slice, 64 x 16 nodes
      std::fprintf(fp, "ix,iz,P_bout\n");
      for (int ix = IX0; ix < IX0 + NXI; ++ix)
        for (int iz = 0; iz < NZ; ++iz)
          std::fprintf(fp, "%d,%d,%.17g\n", ix, iz, PB(ix - IX0, ROW, iz));
      std::fclose(fp);
    }
    size_t nfine = 0;
    if (FILE* fp =
          std::fopen(data_path("transfer_fine_slice.csv").c_str(),
                     "w")) { // the TRANSFERRED field on the 2x mesh, same row
      std::fprintf(fp, "ix_real,iz_real,P_pcms\n");
      for (size_t k = 0; k < tn.size(); ++k)
        if (is_int(tn[k][1]) && std::lround(tn[k][1]) == ROW) {
          std::fprintf(fp, "%.17g,%.17g,%.17g\n", tn[k][0], tn[k][2], tv[k]);
          ++nfine;
        }
      std::fclose(fp);
    }
    std::vector<size_t>
      rc; // 1024 non-coincident fine nodes of that row (iz <= 15)
    for (size_t k = 0; k < tn.size(); ++k)
      if (is_int(tn[k][1]) && std::lround(tn[k][1]) == ROW &&
          tn[k][2] <= NZ - 1 + 1e-12 && !(is_int(tn[k][0]) && is_int(tn[k][2])))
        rc.push_back(k);
    std::mt19937 g2(2024);
    std::shuffle(rc.begin(), rc.end(), g2);
    if (rc.size() > 1024)
      rc.resize(1024);
    Pts<3> tq(rc.size());
    for (size_t c = 0; c < rc.size(); ++c)
      tq[c] = tn[rc[c]];
    const auto fb2 = PB.xz_bilinear(tq);
    if (FILE* fp = std::fopen(data_path("transfer_targets.csv").c_str(), "w")) {
      std::fprintf(fp, "xi_ix,zeta_iz,P_bout,P_pcms,P_eval\n");
      for (size_t c = 0; c < rc.size(); ++c) {
        const size_t k = rc[c];
        std::fprintf(fp, "%.17g,%.17g,%.17g,%.17g,%.17g\n", tn[k][0], tn[k][2],
                     fb2[c], tv[k], ev_src[rc[c]]);
      }
      std::fclose(fp);
    }
    say("PLOT DATA: transfer_slice.csv (%d x %d source nodes), "
        "transfer_fine_slice.csv (%zu fine nodes), transfer_targets.csv (%zu "
        "targets in the iy_file = %d row)\n",
        NXI, NZ, nfine, rc.size(), ROW);
  }

  // Fractional random target nodes (ix, iz)
  // PCMS value versus BOUT++'s XZBilinear
  {
    const int ROW = 32;
    std::mt19937 gr(9001);
    std::uniform_real_distribution<double> rx(IX0, IX0 + NXI - 1 - 1e-9),
      rz(0.0, NZ - 1e-9);
    Pts<3> qr(1024);
    for (auto& p : qr)
      p = {rx(gr), double(ROW), rz(gr)}; // (ix, iy, iz)
    const auto pv = eval_points<3>(src_space, src, qr, IndexSpace3D::Create());
    const auto fb3 = PB.xz_bilinear(qr);
    double d5f = 0;
    if (FILE* fp =
          std::fopen(data_path("transfer_random_targets.csv").c_str(), "w")) {
      std::fprintf(fp, "xi_ix,zeta_iz,P_bout,P_pcms\n");
      for (size_t c = 0; c < qr.size(); ++c) {
        const double b = fb3[c];
        d5f = std::max(d5f, std::fabs(pv[c] - b));
        std::fprintf(fp, "%.17g,%.17g,%.17g,%.17g\n", qr[c][0], qr[c][2], b,
                     pv[c]);
      }
      std::fclose(fp);
    }
    say("[xzbilinear-random] random real (ix, iz) at iy_file = %d, %zu points: "
        "PCMS vs "
        "BOUT++ XZBilinear max diff = %.3e  (relative: %.3e)\n",
        ROW, qr.size(), d5f, d5f / pmax);
  }

  // Nodal RMS and max of the transfer error
  {
    double e2 = 0, emx = 0; // target vs PCMS source evaluation, every node
    for (size_t k = 0; k < tn.size(); ++k) {
      const double e = tv[k] - ev_src[k];
      e2 += e * e;
      emx = std::max(emx, std::fabs(e));
    }
    double c2 = 0, cmx = 0; // target vs BOUT++
    size_t nc = 0;
    for (size_t k = 0; k < tn.size(); ++k)
      if (is_int(tn[k][0]) && is_int(tn[k][1]) && is_int(tn[k][2])) {
        const double e =
          tv[k] - PB(int(std::lround(tn[k][0])) - IX0,
                     int(std::lround(tn[k][1])), int(std::lround(tn[k][2])));
        c2 += e * e;
        cmx = std::max(cmx, std::fabs(e));
        ++nc;
      }
    say("[norm] nodal RMS error of target vs PCMS source evaluation, %zu "
        "nodes: %.3e  (max %.3e)\n",
        tn.size(), std::sqrt(e2 / tn.size()), emx);
    say("[norm] nodal RMS error of target vs BOUT++, %zu coincident nodes: "
        "%.3e  (max %.3e)\n",
        nc, std::sqrt(c2 / nc), cmx);
  }

  say("==== PCMS work inside the pause took %.1f s; returning control to "
      "BOUT++ (the run ends: nout = %d) ====\n",
      secs(t_pause), nout);
}
} // namespace

int main(int argc, char** argv)
{
  Kokkos::ScopeGuard kokkos(argc, argv);
  g_rep = std::fopen(data_path("transfer_result.txt").c_str(), "w");
  bool ran = false;
  const int status = run_elmpb(argc, argv, [] {}, work, ran);
  say("\ntransfer finished: the pause %s, BOUT++ exit status %d\n",
      ran ? "ran" : "did NOT run", status);
  if (g_rep)
    std::fclose(g_rep);
  return ran ? status : 1;
}
