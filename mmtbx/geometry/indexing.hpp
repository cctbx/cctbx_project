#ifndef MMTBX_GEOMETRY_INDEXING_H
#define MMTBX_GEOMETRY_INDEXING_H

#include <scitbx/math/cartesian_product_fixed_size.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/numeric/conversion/converter.hpp>

#include <vector>

namespace mmtbx
{

namespace geometry
{

namespace indexing
{

template< typename Object, typename Algorithm >
class OverlapEqualityFilter : private Algorithm
{
public:
  typedef Object object_type;
  typedef Algorithm algorithm_type;

private:
  object_type object_;

public:
  OverlapEqualityFilter(const Object& object);
  ~OverlapEqualityFilter();

  inline bool operator ()(const object_type& other) const;
};

template< typename Iterator, typename Filter >
struct FilterHelper
{
public:
  typedef boost::filter_iterator< Filter, Iterator > filter_iterator;
  typedef scitbx::math::cartesian_product::iterated_range< filter_iterator >
    filter_range_type;
};

template< typename Iterator, typename PreFilter, typename PostFilter >
struct PrefilterHelper
{
public:
  typedef boost::filter_iterator< PreFilter, Iterator > prefilter_iterator;
  typedef boost::filter_iterator< PostFilter, prefilter_iterator >
    filter_iterator;
  typedef scitbx::math::cartesian_product::iterated_range< filter_iterator >
    filter_range_type;
};

template< typename Object, typename Algorithm >
class Linear
{
public:
  typedef Object object_type;
  typedef Algorithm algorithm_type;
  typedef std::vector< object_type > storage_type;
  typedef OverlapEqualityFilter< object_type, algorithm_type > filter_type;
  typedef FilterHelper< typename storage_type::const_iterator, filter_type >
    filter_helper_type;
  typedef typename filter_helper_type::filter_iterator const_iterator;
  typedef typename filter_helper_type::filter_range_type range_type;

private:
  storage_type objects_;

public:
  Linear();
  ~Linear();

  inline void add(const object_type& object);
  inline range_type overlapping_with(const object_type& object) const;
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

#include "indexing.hxx"

} // namespace indexing
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_INDEXING_H
