template< typename Vector >
Transform< Vector >::Transform(
  const vector_type& centre,
  const value_type& radius
  )
  : centre_( centre ), radius_( radius )
{}

template< typename Vector >
Transform< Vector >::~Transform()
{}

template< typename Vector >
typename Transform< Vector >::vector_type
Transform< Vector>::operator ()(const vector_type& point) const
{
  return radius_ * point + centre_;
}

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

