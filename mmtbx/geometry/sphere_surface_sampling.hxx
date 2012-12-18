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

template< typename Vector >
GoldenSpiral< Vector >::GoldenSpiral(size_t count)
  : count_( count ),
    unit_area_(
      4.0 / count
      * boost::math::constants::pi<
        typename GoldenSpiral< Vector >::value_type
        >()
        )
{
  points_.reserve( count );

  const value_type step = 2.0 / count;
  const value_type offset = step / 2.0;

  for ( size_t index = 0; index < count; ++index)
  {
    const value_type y = step * index - 1.0 + offset;
    const value_type r = std::sqrt( 1.0 - y * y );
    const value_type phi = pitch * index;

    points_.push_back(
      vector_type( std::cos( phi ) * r, y, std::sin( phi ) * r )
      );
  }
}

template< typename Vector >
GoldenSpiral< Vector >::~GoldenSpiral()
{}

template< typename Vector >
typename GoldenSpiral< Vector >::range_type
GoldenSpiral< Vector >::transformed(
  const vector_type& centre,
  const value_type& radius
  ) const
{
  transformer_type t = transformer_type( centre, radius );
  return range_type(
    const_iterator( points_.begin(), t ),
    const_iterator( points_.end(), t )
    );
}

template< typename Vector >
const typename GoldenSpiral< Vector >::value_type&
GoldenSpiral< Vector >::unit_area() const
{
  return unit_area_;
}

template< typename Vector >
const size_t&
GoldenSpiral< Vector >::count() const
{
  return count_;
}

template< typename Vector >
const typename GoldenSpiral< Vector >::storage_type&
GoldenSpiral< Vector >::points() const
{
  return points_;
}

template< typename Vector >
const typename GoldenSpiral< Vector >::value_type
GoldenSpiral< Vector >::pitch =
  boost::math::constants::pi< typename GoldenSpiral< Vector >::value_type >()
  * ( 3.0 - std::sqrt( 5.0 ) );

