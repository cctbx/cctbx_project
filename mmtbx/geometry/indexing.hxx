template< typename Object, typename Vector >
Linear< Object, Vector >::Linear()
{}

template< typename Object, typename Vector >
Linear< Object, Vector >::~Linear()
{}

template< typename Object, typename Vector >
void
Linear< Object, Vector >::add(const object_type& object, const vector_type& position)
{
  objects_.push_back( object );
}

template< typename Object, typename Vector >
typename Linear< Object, Vector >::range_type
Linear< Object, Vector >::close_to(const vector_type& centre) const
{
  return range_type( objects_.begin(), objects_.end() );
}

template< typename Object, typename Vector >
size_t
Linear< Object, Vector >::size() const
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
typename Discretizer< Continuous, Discrete >::continuous_type const&
Discretizer< Continuous, Discrete >::base() const
{
  return base_;
}

template< typename Continuous, typename Discrete >
typename Discretizer< Continuous, Discrete >::discrete_type
Discretizer< Continuous, Discrete >::operator ()(const continuous_type& value)
  const
{
  return converter_type::convert( ( value - base_ ) / unit_ );
}

template< typename Vector, typename Voxel, typename Discrete >
Voxelizer< Vector, Voxel, Discrete >::Voxelizer(
  const vector_type& base,
  const vector_type& step
  )
  : discretizers_(
    discretizer_type( base[0], step[0] ),
    discretizer_type( base[1], step[1] ),
    discretizer_type( base[2], step[2] )
    )
{}

template< typename Vector, typename Voxel, typename Discrete >
Voxelizer< Vector, Voxel, Discrete >::~Voxelizer()
{}

template< typename Vector, typename Voxel, typename Discrete >
typename Voxelizer< Vector, Voxel, Discrete >::voxel_type
Voxelizer< Vector, Voxel, Discrete >::operator ()(const vector_type& vector) const
{
  return voxel_type(
    boost::fusion::at_c<0>( discretizers_ )( vector[0] ),
    boost::fusion::at_c<1>( discretizers_ )( vector[1] ),
    boost::fusion::at_c<2>( discretizers_ )( vector[2] )
    );
}

template< typename Vector, typename Voxel, typename Discrete >
typename Voxelizer< Vector, Voxel, Discrete >::vector_type
Voxelizer< Vector, Voxel, Discrete >::base() const
{
  return vector_type(
    boost::fusion::at_c<0>( discretizers_ ).base(),
    boost::fusion::at_c<1>( discretizers_ ).base(),
    boost::fusion::at_c<2>( discretizers_ ).base()
    );
}

template< typename T >
HashCombine::result_type
HashCombine::operator ()(const T& arg, result_type current) const
{
  boost::hash_combine( current, arg );
  return current;
}

template< typename FusionVector >
typename FusionVectorHasher< FusionVector >::result_type
FusionVectorHasher< FusionVector >::operator ()(const FusionVector& myvec) const
{
  return boost::fusion::fold( myvec, 0, HashCombine() );
}

template< typename Object, typename Vector, typename Discrete >
Hash< Object, Vector, Discrete >::Hash(
  const voxelizer_type& voxelizer,
  const discrete_type& margin
  )
  : voxelizer_( voxelizer ), margin_( margin )
{}

template< typename Object, typename Vector, typename Discrete >
Hash< Object, Vector, Discrete >::~Hash()
{}

template< typename Object, typename Vector, typename Discrete >
void
Hash< Object, Vector, Discrete >::add(const object_type& object, const vector_type& position)
{
  voxel_type key = voxelizer_( position );

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

template< typename Object, typename Vector, typename Discrete >
typename Hash< Object, Vector, Discrete >::range_type
Hash< Object, Vector, Discrete >::close_to(const vector_type& centre) const
{
  range_type result;

  for(
    cartesian_type cit = make_cartesian_iterator_around(
      voxelizer_( centre ),
      voxel_type( 0, 0, 0 )
      );
    cit != cartesian_type::end();
    ++cit
    )
  {
    typename storage_type::const_iterator it = objects_.find( *cit );

    if ( it != objects_.end() )
    {
      result.storage.push_back( bucket_range_type( it->second.begin(), it->second.end() ) );
    }
  }

  return result;
}

template< typename Object, typename Vector, typename Discrete >
typename Hash< Object, Vector, Discrete >::range_type
Hash< Object, Vector, Discrete >::approx_within_sphere(
  const vector_type& centre,
  const distance_type& radius
  ) const
{
  range_type result;

  for(
    cartesian_type cit = make_cartesian_iterator_around(
      voxelizer_( centre ),
      voxelizer_( voxelizer_.base() + vector_type( radius, radius, radius ) )
      );
    cit != cartesian_type::end();
    ++cit
    )
  {
    typename storage_type::const_iterator it = objects_.find( *cit );

    if ( it != objects_.end() )
    {
      result.storage.push_back( bucket_range_type( it->second.begin(), it->second.end() ) );
    }
  }

  return result;
}

template< typename Object, typename Vector, typename Discrete >
size_t
Hash< Object, Vector, Discrete >::size() const
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

template< typename Object, typename Vector, typename Discrete >
size_t
Hash< Object, Vector, Discrete >::cubes() const
{
  return objects_.size();
}

template< typename Object, typename Vector, typename Discrete >
typename Hash< Object, Vector, Discrete >::cartesian_type
Hash< Object, Vector, Discrete >::make_cartesian_iterator_around(
  const voxel_type& voxel,
  const voxel_type& margin
  ) const
{
  typedef typename
    boost::mpl::at_c< typename cartesian_type::range_types_vector, 0 >::type
    cart_range_type;
  const discrete_type& voxel_x = boost::fusion::at_c<0>( voxel );
  const discrete_type& voxel_y = boost::fusion::at_c<1>( voxel );
  const discrete_type& voxel_z = boost::fusion::at_c<2>( voxel );
  const discrete_type& margin_x = boost::fusion::at_c<0>( margin );
  const discrete_type& margin_y = boost::fusion::at_c<1>( margin );
  const discrete_type& margin_z = boost::fusion::at_c<2>( margin );
  typedef typename cartesian_type::range_vector_type initializer_type;
  return cartesian_type(
    initializer_type(
      cart_range_type(
        itercount( voxel_x - margin_x - margin_ ),
        itercount( voxel_x + margin_x + margin_ + 1 )
        ),
      cart_range_type(
        itercount( voxel_y - margin_y - margin_),
        itercount( voxel_y + margin_y + margin_ + 1 )
        ),
      cart_range_type(
        itercount( voxel_z - margin_z - margin_ ),
        itercount( voxel_z + margin_z + margin_ + 1 )
        )
      )
    );
}

