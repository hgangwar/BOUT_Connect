// The interface between the BOUT++ half and the PCMS half of the in-process
// tests: the grid constants, BoutSnapshot (P copied out of BOUT++ memory) and
// the entry points implemented in src/elmpb.cpp. Includes neither BOUT++ nor
// PCMS.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace bout_connect
{

// The elm-pb (cbm18_dens8, nx = 68, ny = 64, nz = 16) index grid shared by
// every test.
constexpr int NZ = 16, NY = 64, NXI = 64, IX0 = 2;

struct BoutSnapshot
{
  std::vector<double>
    pb; // P: pb[iz + NZ * (iy + NY * (ix - IX0))], iy = iy_file
  double pmax = 0, pmin = 1e300, pmaxv = -1e300; // max|P|, min P, max P
  std::string attrs;        // direction_y, direction_z, cell_location
  double dz = 0;            // toroidal step, radians
  std::vector<double> R, Z; // Rxy, Zxy:

  // P at the node (ix - IX0, iy_file, iz); iz wraps periodically.
  double operator()(int ixr, int iy, int iz) const
  {
    return pb[size_t(((iz % NZ) + NZ) % NZ) +
              NZ * (size_t(iy) + size_t(NY) * ixr)];
  }

  // BOUT++'s own XZBilinear
  std::function<std::vector<double>(const std::vector<std::array<double, 3>>&)>
    xz_bilinear;
};

using SnapshotWork =
  std::function<void(const BoutSnapshot&, double t, int iter, int nout)>;

int run_elmpb(int argc, char** argv, const std::function<void()>& after_init,
              const SnapshotWork& work, bool& ran);

int run_elmpb_monitor(
  int argc, char** argv,
  const std::function<void(double t, int iter, int nout)>& on_call);

int option_int(const char* section, const char* key, const char* doc, int dflt);
double option_double(const char* section, const char* key, const char* doc,
                     double dflt);
bool option_bool(const char* section, const char* key, const char* doc,
                 bool dflt);

} // namespace bout_connect
