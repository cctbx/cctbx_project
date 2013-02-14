template< typename Object >
Linear< Object >::Linear()
{}

template< typename Object >
Linear< Object >::~Linear()
{}

template< typename Object >
void
Linear< Object >::add(const object_type& object)
{
  objects_.push_back( object );
}

template< typename Object >
const typename Linear< Object >::range_type&
Linear< Object >::close_to(const object_type& object) const
{
  return objects_;
}

template< typename Object >
size_t
Linear< Object >::size() const
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

template< typename Object, typename Voxelizer >
Hash< Object, Voxelizer >::Hash(const voxelizer_type& voxelizer)
  : voxelizer_( voxelizer )
{}

template< typename Object, typename Voxelizer >
Hash< Object, Voxelizer >::~Hash()
{}

template< typename Object, typename Voxelizer >
void
Hash< Object, Voxelizer >::add(const object_type& object)
{
  voxel_type low = voxelizer_( object.low() );
  voxel_type high = voxelizer_( object.high() );
  assert( low[0] <= high[0] );
  assert( low[1] <= high[1] );
  assert( low[2] <= high[2] );

  for( discrete_type h = low[0]; h < high[0] + 1; ++h )
  {
    for( discrete_type k = low[1]; k < high[1] + 1; ++k )
    {
      for( discrete_type l = low[2]; l < high[2] + 1; ++l )
      {
        voxel_type key = voxel_type( h, k, l );
        typename storage_type::iterator it = objects_.find( key );

        if ( it == objects_.end() )
        {
          std::pair< typename storage_type::iterator, bool > result =
            objects_.insert(
              typename storage_type::value_type( key, bucket_type() )
              );
          assert ( result.second );
          it = result.first;
        }

        it->second.push_back( object );
      }
    }
  }
}

template< typename Object, typename Voxelizer >
typename Hash< Object, Voxelizer >::range_type
Hash< Object, Voxelizer >::close_to(const object_type& object) const
{
  range_type result;

  voxel_type low = voxelizer_( object.low() );
  voxel_type high = voxelizer_( object.high() );
  assert( low[0] <= high[0] );
  assert( low[1] <= high[1] );
  assert( low[2] <= high[2] );

  for( discrete_type h = low[0]; h < high[0] + 1; ++h )
  {
    for( discrete_type k = low[1]; k < high[1] + 1; ++k )
    {
      for( discrete_type l = low[2]; l < high[2] + 1; ++l )
      {
        voxel_type key = voxel_type( h, k, l );
        typename storage_type::const_iterator it = objects_.find( key );

        if ( it == objects_.end() )
        {
          continue;
        }

        result.insert( it->second.begin(), it->second.end() );
      }
    }
  }

  return result;
}

template< typename Object, typename Voxelizer >
size_t
Hash< Object, Voxelizer >::size() const
{
  range_type result;

  for (
    typename storage_type::const_iterator it = objects_.begin();
    it != objects_.end();
    ++it
    )
  {
    result.insert( it->second.begin(), it->second.end() );
  }

  return result.size();
}

template< typename Object, typename Voxelizer >
size_t
Hash< Object, Voxelizer >::cubes() const
{
  return objects_.size();
}

template< typename Object, typename Voxelizer >
size_t
Hash< Object, Voxelizer >::count() const
{
  size_t count = 0;

  for (
    typename storage_type::const_iterator it = objects_.begin();
    it != objects_.end();
    ++it
    )
  {
    count += it->second.size();
  }

  return count;
}

