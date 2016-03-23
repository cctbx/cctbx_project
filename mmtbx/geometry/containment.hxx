// Point in origin-centred shape
template< typename Vector >
bool
PointInOriginSphere::operator ()(
  const Vector& point,
  const typename Vector::value_type& radius_sq
  )
  const
{
  return point.length_sq() < radius_sq;
}

template< typename Vector >
bool
PointInOriginDiamond::operator ()(
  const Vector& point,
  const typename Vector::value_type& radius
  )
  const
{
  typename Vector::value_type sum = fabs( point[0] )
    + fabs( point[1] )
    + fabs( point[2] );
  return sum < radius;
}

template< typename Vector >
bool
PointInOriginCube::operator ()(
  const Vector& point,
  const typename Vector::value_type& radius
  )
  const
{
  return ( fabs( point[0] ) < radius )
    && ( fabs( point[1] ) < radius )
    && ( fabs( point[2] ) < radius );
}

template< typename FloatType >
bool
PointInCube::one_dimensional(
  const FloatType& point,
  const FloatType& low,
  const FloatType& high
  )
  const
{
  return ( low < point ) && ( point < high );
}

template< typename Vector >
bool
PointInCube::operator ()(
  const Vector& point,
  const Vector& low,
  const Vector& high
  )
  const
{
  typedef typename Vector::value_type value_type;

  return one_dimensional( point[0], low[0], high[0] )
    && one_dimensional( point[1], low[1], high[1] )
    && one_dimensional( point[2], low[2], high[2] );
}

// Containment in multiple shapes
template< bool select_overlapped >
template< typename Container >
bool
PurePythagorean< select_overlapped >::operator ()(
  const typename Container::value_type::vector_type& point,
  const Container& neighbours
  ) const
{
  typedef typename Container::value_type::vector_type vector_type;
  typedef typename Container::const_iterator const_iterator;

  for ( const_iterator it = neighbours.begin(); it != neighbours.end(); ++it )
  {
    if ( PointInOriginSphere::operator ()(
        point - it->centre(),
        it->radius_sq()
        )
      )
    {
      return return_value_overlapped;
    }
  }

  return return_value_accessible;
}

template< bool select_overlapped >
template< typename Container >
bool
DiamondPrefilter< select_overlapped >::operator ()(
  const typename Container::value_type::vector_type& point,
  const Container& neighbours
  ) const
{
  typedef typename Container::value_type::vector_type vector_type;
  typedef typename Container::const_iterator const_iterator;

  for ( const_iterator it = neighbours.begin(); it != neighbours.end(); ++it )
  {
    vector_type diff = point - it->centre();

    if ( PointInOriginDiamond::operator ()( diff, it->radius() ) )
    {
      return return_value_overlapped;
    }
  }

  for ( const_iterator it = neighbours.begin(); it != neighbours.end(); ++it )
  {
    vector_type diff = point - it->centre();

    if ( PointInOriginSphere::operator ()( diff, it->radius_sq() ) )
    {
      return return_value_overlapped;
    }
  }

  return return_value_accessible;
}

template< bool select_overlapped >
template< typename Container >
bool
CubePrefilter< select_overlapped >::operator ()(
  const typename Container::value_type::vector_type& point,
  const Container& neighbours
  ) const
{
  typedef typename Container::value_type::vector_type vector_type;
  typedef typename Container::const_iterator const_iterator;

  for ( const_iterator it = neighbours.begin(); it != neighbours.end(); ++it )
  {
    if ( !PointInCube::operator ()( point, it->low(), it->high() ) )
    {
      continue;
    }

    vector_type diff = point - it->centre();

    if ( PointInOriginSphere::operator ()( diff, it->radius_sq() ) )
    {
      return return_value_overlapped;
    }
  }

  return return_value_accessible;
}

// ContainmentFilter
template< typename Neighbour, typename Algorithm >
Checker< Neighbour, Algorithm >::Checker()
{}

template< typename Neighbour, typename Algorithm >
Checker< Neighbour, Algorithm >::~Checker()
{}

template< typename Neighbour, typename Algorithm >
void
Checker< Neighbour, Algorithm >::add(neighbour_type const& object)
{
  neighbours_.push_back( object );
}

template< typename Neighbour, typename Algorithm >
template< typename InputIterator >
void
Checker< Neighbour, Algorithm >::add(InputIterator begin, InputIterator end)
{
  neighbours_.insert( neighbours_.end(), begin, end );
}

template< typename Neighbour, typename Algorithm >
const typename Checker< Neighbour, Algorithm >::storage_type&
Checker< Neighbour, Algorithm >::neighbours() const
{
  return neighbours_;
}

template< typename Neighbour, typename Algorithm >
bool
Checker< Neighbour, Algorithm >::operator ()(const vector_type& point) const
{
  return Algorithm::operator ()( point, neighbours_ );
}

