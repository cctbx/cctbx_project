template< typename Vector >
size_t
Sphere< Vector >::current = 0;

template< typename Vector >
Sphere< Vector >::Sphere(const vector_type& centre, const value_type& radius)
  : primitive::Sphere< Vector >( centre, radius ), index_( current++ )
{}

template< typename Vector >
Sphere< Vector >::~Sphere()
{}

template< typename Vector >
typename Sphere< Vector >::vector_type
Sphere< Vector >::low() const
{
  const value_type& radius = this->radius();
  vector_type radius_vector = vector_type( radius, radius, radius );
  return this->centre() - radius_vector;
}

template< typename Vector >
typename Sphere< Vector >::vector_type
Sphere< Vector >::high() const
{
  const value_type& radius = this->radius();
  vector_type radius_vector = vector_type( radius, radius, radius );
  return this->centre() + radius_vector;
}

template< typename Vector >
const size_t&
Sphere< Vector >::index() const
{
  return index_;
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

template< typename Sphere >
size_t
hash_value(const Sphere& object)
{
  boost::hash< typename Sphere::index_type > hasher;
  return hasher( object.index() );
}

