template< typename Object, typename Algorithm >
OverlapFilter< Object, Algorithm >::OverlapFilter(const object_type& object)
  : object_( object )
{}

template< typename Object, typename Algorithm >
OverlapFilter< Object, Algorithm >::~OverlapFilter()
{}

template< typename Object, typename Algorithm >
bool
OverlapFilter< Object, Algorithm >::operator ()(const object_type& other) const
{
  return Algorithm::operator ()( object_, other );
}

template< typename Object, typename Algorithm >
Linear< Object, Algorithm >::Linear()
{}


template< typename Object, typename Algorithm >
Linear< Object, Algorithm >::~Linear()
{}

template< typename Object, typename Algorithm >
void
Linear< Object, Algorithm >::add(const object_type& object)
{
  objects_.push_back( object );
}

template< typename Object, typename Algorithm >
typename Linear< Object, Algorithm >::range_type
Linear< Object, Algorithm >::overlapping_with(const object_type& object) const
{
  overlap_filter_type filter = overlap_filter_type( object );
  return range_type(
    const_iterator( filter, objects_.begin(), objects_.end() ),
    const_iterator( filter, objects_.end(), objects_.end() )
    );
}

template< typename Object, typename Algorithm >
template< typename PreFilter >
typename PrefilterHelper<
  typename Linear< Object, Algorithm >::storage_type::const_iterator,
  PreFilter,
  typename Linear< Object, Algorithm >::overlap_filter_type
  >::filter_range_type
Linear< Object, Algorithm >::prefiltered_overlapping_with(
  const object_type& object,
  const PreFilter& prefilter
  ) const
{
  typedef PrefilterHelper<
    typename storage_type::const_iterator,
    PreFilter,
    overlap_filter_type
    > prefilter_helper;
  typedef typename prefilter_helper::prefilter_iterator prefilter_iterator;
  typedef typename prefilter_helper::filter_iterator filter_iterator;

  prefilter_iterator begin = prefilter_iterator(
    prefilter,
    objects_.begin(),
    objects_.end()
    );
  prefilter_iterator end = prefilter_iterator(
    prefilter,
    objects_.end(),
    objects_.end()
    );
  overlap_filter_type filter = overlap_filter_type( object );

  return typename prefilter_helper::filter_range_type(
    filter_iterator( filter, begin, end ),
    filter_iterator( filter, end, end )
    );
}

template< typename Object, typename Algorithm >
size_t
Linear< Object, Algorithm >::size() const
{
  return objects_.size();
}

