#ifndef MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H
#define MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H

// The golden spiral is scitbx::math::golden_spiral; this header keeps its
// former name for code outside cctbx that includes it (phasertng). It can go
// once phasertng's copy of calculator.hpp includes scitbx/math/golden_spiral.h.
#include <scitbx/math/golden_spiral.h>

namespace mmtbx
{

namespace geometry
{

namespace sphere_surface_sampling
{

template< typename Vector >
using GoldenSpiral = scitbx::math::golden_spiral< Vector >;

} // namespace sphere_surface_sampling
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H
