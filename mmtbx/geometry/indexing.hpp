#ifndef MMTBX_GEOMETRY_INDEXING_H
#define MMTBX_GEOMETRY_INDEXING_H

#include <scitbx/math/cartesian_product_fixed_size.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <vector>

namespace mmtbx
{

namespace geometry
{

namespace indexing
{

template< typename Object, typename Algorithm >
class OverlapFilter : private Algorithm
{
public:
  typedef Object object_type;
  typedef Algorithm algorithm_type;

private:
  object_type object_;

public:
  OverlapFilter(const Object& object);
  ~OverlapFilter();

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
  typedef OverlapFilter< object_type, algorithm_type > overlap_filter_type;
  typedef FilterHelper<
    typename storage_type::const_iterator,
    overlap_filter_type
    >
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
  template< typename PreFilter >
  inline typename PrefilterHelper<
    typename storage_type::const_iterator,
    PreFilter,
    overlap_filter_type
    >::filter_range_type
    prefiltered_overlapping_with(
      const object_type& object,
      const PreFilter& prefilter
      ) const;
  inline size_t size() const;
};

#include "indexing.hxx"

} // namespace indexing
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_INDEXING_H
