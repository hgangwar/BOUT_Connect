// Index space (ix, iz) or (ix, iy_file, iz) as a PCMS coordinate system.
#pragma once

#include "pcms/field/coordinate_system.hpp"

#include <memory>
#include <string_view>

namespace bout_connect
{

template <int D>
class IndexSpace final : public pcms::CoordinateSystem
{
  static_assert(D == 2 || D == 3, "index space is 2D (ix, iz) or 3D");

public:
  static std::shared_ptr<const pcms::CoordinateSystem> Create()
  {
    return std::shared_ptr<const pcms::CoordinateSystem>(new IndexSpace());
  }
  std::string_view Kind() const noexcept override
  {
    return D == 3 ? "index_space (ix, iy_file, iz)" : "index_space (ix, iz)";
  }
  int Dimension() const noexcept override { return D; }
  bool HasOrthogonalBasis() const noexcept override { return false; }
  bool HasUnitScaleFactors() const noexcept override { return false; }
  bool operator==(const pcms::CoordinateSystem& o) const noexcept override
  {
    return dynamic_cast<const IndexSpace*>(&o) != nullptr;
  }

private:
  IndexSpace() noexcept = default;
};

using IndexSpace2D = IndexSpace<2>;
using IndexSpace3D = IndexSpace<3>;

} // namespace bout_connect
