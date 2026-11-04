
#pragma once

#include "bout_connect/bout_interface.hpp"
#include "bout_connect/index_space.hpp"

#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/uniform_grid.h"
#include <Kokkos_Core.hpp>

#include <array>
#include <memory>
#include <type_traits>
#include <vector>

namespace bout_connect
{

template <int D>
using Pts = std::vector<std::array<double, D>>;
using CsysPtr = std::shared_ptr<const pcms::CoordinateSystem>;

// Rank-2 PCMS view (device or host) -> host points.
template <int D, class MD>
Pts<D> host_pts(const MD& d)
{
  using VL = std::conditional_t<
    std::is_same_v<typename MD::layout_type, Kokkos::layout_left>,
    Kokkos::LayoutLeft, Kokkos::LayoutRight>;
  Kokkos::View<const double**, VL, pcms::DeviceMemorySpace,
               Kokkos::MemoryUnmanaged>
    dv(d.data_handle(), d.extent(0), d.extent(1));
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dv);
  Pts<D> p(h.extent(0));
  for (size_t k = 0; k < p.size(); ++k)
    for (int i = 0; i < D; ++i)
      p[k][i] = h(k, i);
  return p;
}

// DOF-holder (node) coordinates of a function space.  DOF order is never
// assumed: fill and read values through these coordinates.
template <int D, class S>
Pts<D> node_coords(const S& sp)
{
  return host_pts<D>(sp.GetLayout()->GetDOFHolderCoordinates().GetValues());
}

template <class F>
void set_dofs(F& fn, const std::vector<double>& v)
{
  fn.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      v.data(), pcms::LO(v.size()), 1));
}

template <class F>
std::vector<double> get_dofs(const F& fn, size_t n)
{
  const auto vv = fn.GetDOFHolderDataHost().GetValues();
  std::vector<double> o(n);
  for (size_t k = 0; k < n; ++k)
    o[k] = vv(k, 0);
  return o;
}

// Scalar field evaluated at the points q, given in coordinate system cs.
template <int D, class F>
std::vector<double> eval_points(
  const std::shared_ptr<pcms::LagrangeFunctionSpace>& sp, F& fn,
  const Pts<D>& q, const CsysPtr& cs)
{
  const size_t n = q.size();
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> qq("pts", n, D);
  pcms::CopyHostRank2ViewToDeviceView(
    qq, pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
          q[0].data(), pcms::LO(n), D));
  pcms::CoordinateView<pcms::DeviceMemorySpace> cv(
    cs, pcms::MakeConstRank2View(qq));
  auto ev = sp->CreatePointEvaluator<pcms::Real>(
    pcms::EvaluationRequest::FromCoordinates(cv));
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("out", n, 1);
  ev->Evaluate(fn, pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                     fn.GetData().GetValueBasis(), pcms::MakeRank2View(out)));
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);
  std::vector<double> r(n);
  for (size_t k = 0; k < n; ++k)
    r[k] = h(k, 0);
  return r;
}

// The 3D index-space grid, r times finer than the elm-pb nodes: axis order
// (ix, iy, iz), iz plane NZ is the periodic copy of plane 0.  The first two
// divisions are EQUAL on purpose (63, 63): PCMS's 3D index helper
// (GetDimensionedIndex) divides the stride by divisions[0] where
// divisions[1] is needed, so unequal first two divisions return wrong values.
inline std::shared_ptr<pcms::LagrangeFunctionSpace> make_index_space(int r = 1)
{
  pcms::UniformGrid<3> g;
  g.bot_left = {double(IX0), 0.0, 0.0};
  g.edge_length = {double(NXI - 1), double(NY - 1), double(NZ)};
  g.divisions = {(NXI - 1) * r, (NY - 1) * r, NZ * r};
  return pcms::LagrangeFunctionSpace::FromUniformGrid(
    g, 1, IndexSpace3D::Create(), 1);
}

// The ix-iz slice (ix 2..65, iz 0..16) at fixed iy: 2D grid with UNEQUAL
// divisions (63 x 16 cells, times r), in coordinate system cs.
inline std::shared_ptr<pcms::LagrangeFunctionSpace> make_slice_space(
  const CsysPtr& cs, int r = 1)
{
  pcms::UniformGrid<2> g;
  g.bot_left = {double(IX0), 0.0};
  g.edge_length = {double(NXI - 1), double(NZ)};
  g.divisions = {(NXI - 1) * r, NZ * r};
  return pcms::LagrangeFunctionSpace::FromUniformGrid(g, 1, cs, 1);
}

} // namespace bout_connect
