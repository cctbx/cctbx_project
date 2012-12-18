#ifndef MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H
#define MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H

#include <scitbx/math/cartesian_product_fixed_size.hpp>

#include <boost/iterator/transform_iterator.hpp>
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
class Transform
{
public:
  typedef Vector vector_type;
  typedef typename Vector::value_type value_type;
  typedef vector_type result_type;

private:
  vector_type centre_;
  value_type radius_;

public:
  Transform(const vector_type& centre, const value_type& radius);
  ~Transform();

  inline vector_type operator ()(const vector_type& point) const;
};

template< typename Vector >
class GoldenSpiral
{
public:
  typedef Vector vector_type;
  typedef typename Vector::value_type value_type;
  typedef std::vector< vector_type > storage_type;
  typedef typename storage_type::const_iterator const_storage_iterator;
  typedef Transform< Vector > transformer_type;
  typedef boost::transform_iterator< transformer_type, const_storage_iterator >
    const_iterator;
  typedef scitbx::math::cartesian_product::iterated_range< const_iterator >
    range_type;

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
  inline range_type transformed(
    const vector_type& centre,
    const value_type& radius
    )
    const;
};

#include "sphere_surface_sampling.hxx"

} // namespace sphere_surface_sampling
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_SPHERE_SURFACE_SAMPLING_H

