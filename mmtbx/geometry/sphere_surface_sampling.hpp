#ifndef MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H
#define MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H

#include <boost/math/constants/constants.hpp>

#include <vector>
#include <cmath>

namespace mmtbx
{

namespace geometry
{

namespace sphere_surface_sampling
{

template< typename Vector >
class GoldenSpiral
{
public:
  typedef Vector vector_type;
  typedef typename Vector::value_type value_type;
  typedef std::vector< vector_type > storage_type;
  typedef typename storage_type::const_iterator const_storage_iterator;

private:
  size_t count_;
  value_type unit_area_;
  storage_type points_;

  static const value_type pitch;

public:
  GoldenSpiral(size_t count);
  ~GoldenSpiral();

  inline const value_type& unit_area() const;
  inline const size_t& count() const;
  inline const storage_type& points() const;
};

#include "sphere_surface_sampling.hxx"

} // namespace sphere_surface_sampling
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H

