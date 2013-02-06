template< typename Object, typename Algorithm >
OverlapEqualityFilter< Object, Algorithm >::OverlapEqualityFilter(
  const object_type& object
  )
  : object_( object )
{}

template< typename Object, typename Algorithm >
OverlapEqualityFilter< Object, Algorithm >::~OverlapEqualityFilter()
{}

template< typename Object, typename Algorithm >
bool
OverlapEqualityFilter< Object, Algorithm >::operator ()(
  const object_type& other
  )
  const
{
  return ( other != object_ ) && Algorithm::operator ()( object_, other );
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
  filter_type filter = filter_type( object );
  return range_type(
    const_iterator( filter, objects_.begin(), objects_.end() ),
    const_iterator( filter, objects_.end(), objects_.end() )
    );
}

template< typename Object, typename Algorithm >
size_t
Linear< Object, Algorithm >::size() const
{
  return objects_.size();
}

template< typename Continuous, typename Discrete >
Discretizer< Continuous, Discrete >::Discretizer(
  const continuous_type& base,
  const continuous_type& unit
  )
  : base_( base ), unit_( unit )
{
  assert ( 0 < unit_ );
}

template< typename Continuous, typename Discrete >
Discretizer< Continuous, Discrete >::~Discretizer()
{}

template< typename Continuous, typename Discrete >
typename Discretizer< Continuous, Discrete >::discrete_type
Discretizer< Continuous, Discrete >::operator ()(const continuous_type& value)
  const
{
  return converter_type::convert( ( value - base_ ) / unit_ );
}

template< typename Vector, typename Voxel >
Voxelizer< Vector, Voxel >::Voxelizer(
  const vector_type& base,
  const vector_type& step
  )
  : discretizers_(
    discretizer_type( base[0], step[0] ),
    discretizer_type( base[1], step[1] ),
    discretizer_type( base[2], step[2] )
    )
{}

template< typename Vector, typename Voxel >
Voxelizer< Vector, Voxel >::~Voxelizer()
{}

template< typename Vector, typename Voxel >
typename Voxelizer< Vector, Voxel >::voxel_type
Voxelizer< Vector, Voxel >::operator ()(const vector_type& vector) const
{
  return voxel_type(
    boost::fusion::at_c<0>( discretizers_ )( vector[0] ),
    boost::fusion::at_c<1>( discretizers_ )( vector[1] ),
    boost::fusion::at_c<2>( discretizers_ )( vector[2] )
    );
}

