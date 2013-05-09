#ifndef MMTBX_GEOMETRY_OVERLAP_H
#define MMTBX_GEOMETRY_OVERLAP_H

#include <cmath>

namespace mmtbx
{

namespace geometry
{

namespace overlap
{

struct BetweenSpheres
{
  template< typename Sphere >
  inline bool operator ()(const Sphere& left, const Sphere& right) const;
};

struct BetweenSpheresTolerance
{
  template< typename Sphere >
  inline bool operator ()(
    const Sphere& left,
    const Sphere& right,
    const typename Sphere::value_type& tolerance
    ) const;
};

struct BetweenBoxes
{
  template< typename FloatType >
  inline bool one_dimensional_overlap(
    const FloatType& l_low,
    const FloatType& l_high,
    const FloatType& r_low,
    const FloatType& r_high
    ) const ;
  template< typename Box >
  inline bool operator ()(const Box& left, const Box& right) const;
};

#include "overlap.hxx"

} // namespace overlap
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_OVERLAP_H

