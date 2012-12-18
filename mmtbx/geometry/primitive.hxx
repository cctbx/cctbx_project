// Sphere
template< typename Vector >
Sphere< Vector >::Sphere(const vector_type& centre, const value_type& radius)
  : centre_( centre ), radius_( radius ), radius_sq_( radius * radius )
{}

template< typename Vector >
Sphere< Vector >::~Sphere()
{}

template< typename Vector >
const typename Sphere< Vector >::vector_type&
Sphere< Vector >::centre() const
{
  return centre_;
}

template< typename Vector >
const typename Sphere< Vector >::value_type&
Sphere< Vector >::radius() const
{
  return radius_;
}

template< typename Vector >
const typename Sphere< Vector >::value_type&
Sphere< Vector >::radius_sq() const
{
  return radius_sq_;
}

// Box
template< typename Vector >
Box< Vector >::Box(const vector_type& low, const vector_type& high)
  : low_( low ), high_( high )
{}

template< typename Vector >
Box< Vector >::~Box()
{}

template< typename Vector >
const typename Box< Vector >::vector_type&
Box< Vector >::low() const
{
  return low_;
}

template< typename Vector >
const typename Box< Vector >::vector_type&
Box< Vector >::high() const
{
  return high_;
}

template< typename Vector >
Box< Vector >
Box< Vector >::from_corners(
  const vector_type& corner1,
  const vector_type& corner2
  )
{
  vector_type low =  corner1;
  vector_type high = corner2;
  low.each_update_min( corner2 );
  high.each_update_max( corner1 );

  return Box( low, high );
}

template< typename Vector >
Box< Vector >
Box< Vector >::around_sphere(
  const vector_type& centre,
  const value_type& radius
  )
{
  vector_type radius_vector = vector_type( radius, radius, radius );

  return Box( centre - radius_vector, centre + radius_vector );
}

template< typename Vector >
Box< Vector >
Box< Vector >::from_sphere(const Sphere< Vector >& sphere)
{
  const vector_type& centre = sphere.centre();
  const value_type& radius = sphere.radius();
  vector_type radius_vector = vector_type( radius, radius, radius );

  return Box( centre - radius_vector, centre + radius_vector );
}

// BSphere
template< typename Vector >
BSphere< Vector >::BSphere(const vector_type& centre, const value_type& radius)
  : Sphere< Vector >( centre, radius ),
    Box< Vector >( Box< Vector >::around_sphere( centre, radius ) )
{}

template< typename Vector >
BSphere< Vector >::~BSphere()
{}

