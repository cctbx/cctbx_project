template< typename Vector >
size_t
Sphere< Vector >::current = 0;

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
Sphere< Vector >
Sphere< Vector >::create(const vector_type& centre, const value_type& radius)
{
  return Sphere< Vector >( centre, radius, current++ );
}

template< typename Vector >
bool
operator ==(const Sphere< Vector >& left, const Sphere< Vector >& right)
{
  return left.index() == right.index();
}

template< typename Vector >
bool
operator !=(const Sphere< Vector >& left, const Sphere< Vector >& right)
{
  return left.index() != right.index();
}

