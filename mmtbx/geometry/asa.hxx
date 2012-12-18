template< typename Vector >
Sphere< Vector >::Sphere(
  const vector_type& centre,
  const value_type& radius,
  size_t index
  )
  : primitive::Sphere< Vector >( centre, radius ), index_( index )
{}

template< typename Vector >
Sphere< Vector >::~Sphere()
{}

template< typename Vector >
const size_t&
Sphere< Vector >::index() const
{
  return index_;
}

template< typename Vector >
IndexFilter< Vector >::IndexFilter(const size_t& index)
  : index_( index )
{}

template< typename Vector >
IndexFilter< Vector >::~IndexFilter()
{}

template< typename Vector >
bool
IndexFilter< Vector >::operator ()(const Sphere< Vector >& object) const
{
  return object.index() != index_;
}

