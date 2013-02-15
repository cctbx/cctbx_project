#ifndef MMTBX_GEOMETRY_INDEXING_H
#define MMTBX_GEOMETRY_INDEXING_H

#include <scitbx/math/cartesian_product_fixed_size.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <boost/numeric/conversion/converter.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>

#include <vector>

namespace mmtbx
{

namespace geometry
{

namespace indexing
{

template< typename Object >
class Linear
{
public:
  typedef Object object_type;
  typedef std::vector< object_type > storage_type;
  typedef storage_type range_type;

private:
  storage_type objects_;

public:
  Linear();
  ~Linear();

  inline void add(const object_type& object);
  inline const range_type& close_to(const object_type& object) const;
  inline size_t size() const;
};

template< typename Continuous, typename Discrete >
class Discretizer
{
public:
  typedef Continuous continuous_type;
  typedef Discrete discrete_type;
  typedef boost::numeric::converter<
    Discrete,
    Continuous,
    boost::numeric::conversion_traits< Discrete, Continuous >,
    boost::numeric::def_overflow_handler,
    boost::numeric::Floor< Continuous >
    > converter_type;

private:
  continuous_type base_;
  continuous_type unit_;

public:
  Discretizer(const continuous_type& base, const continuous_type& unit);
  ~Discretizer();

  discrete_type operator ()(const continuous_type& value) const;
};

template< typename Vector, typename Voxel >
class Voxelizer
{
public:
  typedef Vector vector_type;
  typedef Voxel voxel_type;
  typedef typename vector_type::value_type continuous_type;
  typedef typename voxel_type::value_type discrete_type;
  typedef Discretizer< continuous_type, discrete_type > discretizer_type;
  typedef boost::fusion::vector3<
    discretizer_type, discretizer_type, discretizer_type
    > discretizer_vector_type;

private:
  discretizer_vector_type discretizers_;

public:
  Voxelizer(const vector_type& base, const vector_type& step);
  ~Voxelizer();

  voxel_type operator ()(const vector_type& vector) const;
};

template< typename Object, typename Voxelizer >
class Hash
{
public:
  typedef Object object_type;
  typedef Voxelizer voxelizer_type;
  typedef typename object_type::vector_type vector_type;
  typedef typename voxelizer_type::voxel_type voxel_type;
  typedef typename voxelizer_type::discrete_type discrete_type;
  typedef std::vector< object_type > bucket_type;
  typedef boost::unordered_map< voxel_type, bucket_type > storage_type;
  typedef boost::unordered_set< object_type > range_type;

private:
  typedef boost::counting_iterator< discrete_type > itercount;
  typedef boost::mpl::vector< itercount, itercount, itercount > itercount_list;
  typedef scitbx::math::cartesian_product::fixed_size_iterator< itercount_list >
    cartesian_type;

private:
  voxelizer_type voxelizer_;
  storage_type objects_;

public:
  Hash(const voxelizer_type& voxelizer);
  ~Hash();

  inline void add(const object_type& object);
  inline range_type close_to(const object_type& object) const;
  inline size_t size() const;
  inline size_t cubes() const;
  inline size_t count() const;

private:
  cartesian_type make_cartesian_iterator(
    const voxel_type& low,
    const voxel_type& high
    ) const;
};

template< typename OutputVector >
struct vector_reformat
{
  template< typename InputVector >
  OutputVector operator ()(const InputVector& input);
};

#include "indexing.hxx"

} // namespace indexing
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_INDEXING_H
