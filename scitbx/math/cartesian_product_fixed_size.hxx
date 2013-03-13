template< class Iterator >
counter< Iterator >::counter()
  : current( range.begin() )
{}

template< class Iterator >
counter< Iterator >::counter(const range_type& range_)
  : range( range_ ), current( range.begin() )
{}

template< class Iterator >
counter< Iterator >::~counter()
{}

template< class Iterator >
bool
operator ==(const counter< Iterator >& left, const counter< Iterator >& right)
{
  return ( left.current == right.current ) && ( left.range == right.range );
}

template< typename CounterVector, unsigned Depth >
bool
increment_fast_back< CounterVector, Depth >::process(CounterVector& cv) const
{
  typedef typename boost::remove_reference<
    typename boost::fusion::result_of::at_c< CounterVector, Depth >::type
    >::type counter_type;
  counter_type& mycounter = boost::fusion::at_c< Depth >( cv );
  ++( mycounter.current );

  if ( mycounter.current == mycounter.range.end() )
  {
    mycounter.current = mycounter.range.begin();
    return increment_fast_back< CounterVector, Depth - 1 >().process( cv );
  }
  else
  {
    return false;
  }
}

template< typename CounterVector >
bool
increment_fast_back< CounterVector, 0 >::process(CounterVector& cv) const
{
  typedef typename boost::remove_reference<
    typename boost::fusion::result_of::at_c< CounterVector, 0 >::type
    >::type counter_type;
  counter_type& mycounter = boost::fusion::at_c< 0 >( cv );
  ++( mycounter.current );

  return mycounter.current == mycounter.range.end();
}

template< typename CounterVector, unsigned Depth >
bool
increment_fast_front< CounterVector, Depth >::process(CounterVector& cv) const
{
  typedef typename boost::fusion::result_of::size< CounterVector >::type
    mpl_integral_constant;
  typedef typename boost::remove_reference<
    typename boost::fusion::result_of::at_c<
      CounterVector,
      mpl_integral_constant::value - Depth - 1
      >::type
    >::type counter_type;
  counter_type& mycounter = boost::fusion::at_c<
    mpl_integral_constant::value - Depth - 1
    >( cv );
  ++( mycounter.current );

  if ( mycounter.current == mycounter.range.end() )
  {
    mycounter.current = mycounter.range.begin();
    return increment_fast_front< CounterVector, Depth - 1 >().process( cv );
  }
  else
  {
    return false;
  }
}

template< typename CounterVector >
bool
increment_fast_front< CounterVector, 0 >::process(CounterVector& cv) const
{
  typedef typename boost::fusion::result_of::size< CounterVector >::type
    mpl_integral_constant;
  typedef typename boost::remove_reference<
    typename boost::fusion::result_of::at_c<
      CounterVector,
      mpl_integral_constant::value - 1
      >::type
    >::type counter_type;
  counter_type& mycounter = boost::fusion::at_c<
    mpl_integral_constant::value - 1
    >( cv );
  ++( mycounter.current );

  return mycounter.current == mycounter.range.end();
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
fixed_size_iterator< IteratorTypeList, IterationOrder >::fixed_size_iterator()
  : finished_( true )
{}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
fixed_size_iterator< IteratorTypeList, IterationOrder >::fixed_size_iterator(
  const range_vector_type& rv
  )
  : counter_vector_( boost::fusion::transform( rv, get_counter_for_range() ) ),
    finished_( boost::fusion::any( counter_vector_, is_counter_empty() ) )
{
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
fixed_size_iterator< IteratorTypeList, IterationOrder >::~fixed_size_iterator()
{}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
typename fixed_size_iterator< IteratorTypeList, IterationOrder >::iterator_type&
fixed_size_iterator< IteratorTypeList, IterationOrder >::operator ++()
{
  finished_ = order_type::process( counter_vector_ );
  return *this;
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
typename fixed_size_iterator< IteratorTypeList, IterationOrder >::iterator_type
fixed_size_iterator< IteratorTypeList, IterationOrder >::operator ++(int)
{
  iterator_type tmp( *this );
  operator ++();
  return tmp;
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
bool
fixed_size_iterator< IteratorTypeList, IterationOrder >::operator ==(
  const iterator_type& other
  ) const
{
  if ( finished_  || other.finished_ )
  {
    return finished_ == other.finished_;
  }
  else
  {
    return counter_vector_ == other.counter_vector_;
  }
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
bool
fixed_size_iterator< IteratorTypeList, IterationOrder >::operator !=(
  const iterator_type& other
  ) const
{
  return !( *this == other );
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
typename fixed_size_iterator< IteratorTypeList, IterationOrder >::value_type&
fixed_size_iterator< IteratorTypeList, IterationOrder >::operator *()
{
    current_value_ = boost::fusion::transform(
      counter_vector_,
      get_counter_current_value()
      );

    return current_value_;
}

template<
  typename IteratorTypeList,
  template< class, unsigned > class IterationOrder
  >
fixed_size_iterator< IteratorTypeList, IterationOrder >
fixed_size_iterator< IteratorTypeList, IterationOrder >::end()
{
  return iterator_type();
}

