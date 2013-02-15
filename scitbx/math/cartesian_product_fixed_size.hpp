#ifndef SCITBX_MATH_CARTESIAN_PRODUCT_FIXED_H
#define SCITBX_MATH_CARTESIAN_PRODUCT_FIXED_H

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/include/transform.hpp>
#include <boost/fusion/sequence/comparison/equal_to.hpp>
#include <boost/fusion/include/equal_to.hpp>
#include <boost/fusion/algorithm/query/any.hpp>
#include <boost/fusion/include/any.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/include/size.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/size.hpp>

#include <boost/type_traits/remove_reference.hpp>
#include <boost/range/iterator_range.hpp>

namespace scitbx
{

namespace math
{

namespace cartesian_product
{

template< class Iterator >
struct counter
{
  typedef Iterator iterator_type;
  typedef typename iterator_type::value_type value_type;
  typedef boost::iterator_range< iterator_type > range_type;
  typedef counter< iterator_type > counter_type;

  range_type range;
  iterator_type current;

  counter();
  counter(const range_type& range);
  ~counter();
};

template< class Iterator >
bool operator ==(
  const counter< Iterator >& left,
  const counter< Iterator >& right
  );

template< typename CounterVector, unsigned Depth >
struct increment_fast_back
{
  bool process(CounterVector& cv) const;
};

template< typename CounterVector >
struct increment_fast_back< CounterVector, 0 >
{
  bool process(CounterVector& cv) const;
};

template< typename CounterVector, unsigned Depth >
struct increment_fast_front
{
  bool process(CounterVector& cv) const;
};

template< typename CounterVector >
struct increment_fast_front< CounterVector, 0 >
{
  bool process(CounterVector& cv) const;
};

template< typename IteratorTypeList >
class fixed_size_iterator_helper
{
private:
  struct get_value_type
  {
    template< class Any >
    struct apply
    {
      typedef typename Any::value_type type;
    };
  };

  struct get_counter_type
  {
    template< class Iterator >
    struct apply
    {
      typedef typename counter< Iterator >::counter_type type;
    };
  };

  struct get_range_type
  {
    template< class Counter >
    struct apply
    {
      typedef typename Counter::range_type type;
    };
  };

public:
  typedef typename boost::mpl::transform<
    IteratorTypeList,
    get_value_type
    >::type mpl_value_types_vector;
  typedef typename boost::mpl::transform<
    IteratorTypeList,
    get_counter_type
    >::type mpl_counter_types_vector;
  typedef typename boost::mpl::transform<
    mpl_counter_types_vector,
    get_range_type
    >::type mpl_range_types_vector;

public:
  typedef typename boost::fusion::result_of::as_vector<
    mpl_value_types_vector
    >::type value_type;
  typedef typename boost::fusion::result_of::as_vector<
    mpl_counter_types_vector
    >::type counter_vector_type;
  typedef typename boost::fusion::result_of::as_vector<
    mpl_range_types_vector
    >::type range_vector_type;
};

template<
  typename IteratorTypeList,
  template<class, unsigned> class IterationOrder = increment_fast_back
  >
class fixed_size_iterator : IterationOrder<
  typename fixed_size_iterator_helper< IteratorTypeList >::counter_vector_type,
  boost::mpl::size< IteratorTypeList >::value - 1
  >
{
private:
  struct is_counter_empty
  {
    template<typename Counter>
    bool operator()(const Counter& counter) const
    {
      return counter.range.begin() == counter.range.end();
    }
  };

  struct get_counter_for_range
  {
    template< typename Signature > struct result;

    template< typename Transformation, class Range >
    struct result< Transformation( Range ) >
    {
      typedef typename Range::iterator iterator_type;
      typedef typename counter< iterator_type >::counter_type type;
    };

    template< typename Transformation, class Range >
    struct result< Transformation( Range& ) >
    {
      typedef typename Range::iterator iterator_type;
      typedef typename counter< iterator_type >::counter_type type;
    };

    template< class Range >
    typename result< get_counter_for_range( Range ) >::type
    operator()(const Range& range) const
    {
      return typename result< get_counter_for_range( Range ) >::type( range );
    }
  };

  struct get_counter_current_value
  {
    template< typename Signature > struct result;

    template< typename Transformation, class Counter >
    struct result< Transformation( Counter ) >
    {
      typedef typename Counter::value_type type;
    };

    template< typename Transformation, class Counter >
    struct result< Transformation( Counter& ) >
    {
      typedef typename Counter::value_type type;
    };

    template< class Counter >
    typename result< get_counter_current_value( Counter ) >::type
    operator()(const Counter& counter) const
    {
      return *counter.current;
    }
  };

public:
  typedef fixed_size_iterator< IteratorTypeList, IterationOrder > iterator_type;
  typedef fixed_size_iterator_helper< IteratorTypeList > types_helper;
  typedef typename types_helper::mpl_value_types_vector value_types_vector;
  typedef typename types_helper::mpl_counter_types_vector counter_types_vector;
  typedef typename types_helper::mpl_range_types_vector range_types_vector;

  typedef typename types_helper::value_type value_type;
  typedef typename types_helper::counter_vector_type counter_vector_type;
  typedef typename types_helper::range_vector_type range_vector_type;

  typedef IterationOrder<
    counter_vector_type,
    boost::mpl::size< IteratorTypeList >::value - 1
    > order_type;

private:
  counter_vector_type counter_vector_;
  bool finished_;
  value_type current_value_;

public:
  fixed_size_iterator();
  fixed_size_iterator(const range_vector_type& rv);
  ~fixed_size_iterator();

  iterator_type& operator ++();
  iterator_type operator ++(int);
  bool operator ==(const iterator_type& rhs) const;
  bool operator !=(const iterator_type& rhs) const;
  value_type& operator *();

  static iterator_type end();
};

template<
  typename ContainerTypeList,
  template<class, unsigned> class IterationOrder = increment_fast_back
  >
class from_containers
{
private:
  struct get_iterator_type
  {
    template< class Container >
    struct apply
    {
      typedef typename Container::const_iterator type;
    };
  };

public:
  typedef typename boost::mpl::transform<
    ContainerTypeList,
    get_iterator_type
    >::type mpl_iterator_types_vector;
  typedef fixed_size_iterator< mpl_iterator_types_vector, IterationOrder >
    cartesian_iterator_type;
  typedef typename cartesian_iterator_type::range_vector_type
    initializer_type;
};

struct get_range_from_container
{
  template< typename Signature > struct result;

  template< typename Transformation, class Container >
  struct result< Transformation( Container ) >
  {
    typedef typename Container::const_iterator iterator_type;
    typedef typename counter< iterator_type >::range_type type;
  };

  template< class Container >
  typename result< get_range_from_container( Container ) >::type
  operator()(const Container& container) const
  {
    typedef typename result< get_range_from_container( Container ) >::type
      result_type;
    return result_type( container.begin(), container.end() );
  }
};

#include "cartesian_product_fixed_size.hxx"

} // namespace cartesian_product
} // namespace geometry
} // namespace mmtbx

#endif // SCITBX_MATH_CARTESIAN_PRODUCT_FIXED_H

