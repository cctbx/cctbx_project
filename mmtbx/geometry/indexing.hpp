#ifndef MMTBX_GEOMETRY_INDEXING_H
#define MMTBX_GEOMETRY_INDEXING_H

#include <mmtbx/geometry/flattening.hpp>
#include <scitbx/math/cartesian_product_fixed_size.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <boost/numeric/conversion/converter.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>

#include <boost/mpl/vector.hpp>

#include <vector>

namespace mmtbx
{

namespace geometry
{

namespace indexing
{

template< typename Object, typename Vector >
class Linear
{
public:
  typedef Object object_type;
  typedef Vector vector_type;
  typedef std::vector< object_type > storage_type;
  typedef boost::iterator_range< typename storage_type::const_iterator > range_type;

private:
  storage_type objects_;

public:
  Linear();
  ~Linear();

  inline void add(const object_type& object, const vector_type& position);
  inline range_type close_to(const vector_type& centre) const;
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

  continuous_type const& base() const;
  discrete_type operator ()(const continuous_type& value) const;
};

template< typename Vector, typename Voxel, typename Discrete >
class Voxelizer
{
public:
  typedef Vector vector_type;
  typedef Voxel voxel_type;
  typedef typename vector_type::value_type continuous_type;
  typedef Discrete discrete_type;
  typedef Discretizer< continuous_type, discrete_type > discretizer_type;
  typedef boost::fusion::vector3<
    discretizer_type, discretizer_type, discretizer_type
    > discretizer_vector_type;

private:
  discretizer_vector_type discretizers_;

public:
  Voxelizer(const vector_type& base, const vector_type& step);
  ~Voxelizer();

  vector_type base() const;
  voxel_type operator ()(const vector_type& vector) const;
};

struct HashCombine
{
  typedef std::size_t result_type;

  template< typename T >
  result_type operator ()(const T& arg, result_type current) const;
};

template< typename FusionVector >
struct FusionVectorHasher
{
  typedef std::size_t result_type;
  result_type operator ()(const FusionVector& myvector) const;
};

template< typename Object, typename Vector, typename Discrete >
class Hash
{
public:
  typedef Object object_type;
  typedef Vector vector_type;
  typedef Discrete discrete_type;

private:
  typedef boost::counting_iterator< discrete_type > itercount;
  typedef boost::mpl::vector< itercount, itercount, itercount > itercount_list;
  typedef scitbx::math::cartesian_product::fixed_size_iterator< itercount_list >
    cartesian_type;

public:
  typedef typename vector_type::value_type distance_type;
  typedef typename cartesian_type::value_type voxel_type;
  typedef Voxelizer< vector_type, voxel_type, discrete_type > voxelizer_type;
  typedef std::vector< object_type > bucket_type;
  typedef FusionVectorHasher< voxel_type > hasher_type;
  typedef boost::unordered_map< voxel_type, bucket_type, hasher_type > storage_type;
  typedef boost::iterator_range< typename bucket_type::const_iterator >
    bucket_range_type;
  typedef utility::flattening_range< bucket_range_type > range_type;

private:
  voxelizer_type voxelizer_;
  storage_type objects_;
  discrete_type margin_;

public:
  Hash(const voxelizer_type& voxelizer, const discrete_type& margin);
  ~Hash();

  inline void add(const object_type& object, const vector_type& position);
  inline range_type close_to(const vector_type& centre) const;
  inline range_type approx_within_sphere(
    const vector_type& centre,
    const distance_type& radius
    ) const;
  inline size_t size() const;
  inline size_t cubes() const;

private:
  cartesian_type make_cartesian_iterator_around(
    const voxel_type& voxel,
    const voxel_type& margin
    ) const;
};

#include "indexing.hxx"

} // namespace indexing
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_INDEXING_H
