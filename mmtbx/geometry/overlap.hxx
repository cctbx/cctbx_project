template< typename Sphere >
bool
BetweenSpheres::operator ()(const Sphere& left, const Sphere& right) const
{
  typedef typename Sphere::value_type value_type;
  value_type sum_radii = left.radius() + right.radius();
  return ( left.centre() - right.centre() ).length_sq() < sum_radii * sum_radii;
}

template< typename FloatType >
bool
BetweenBoxes::one_dimensional_overlap(
  const FloatType& l_low,
  const FloatType& l_high,
  const FloatType& r_low,
  const FloatType& r_high
  ) const
{
  FloatType distance = fabs( l_low + l_high - r_low - r_high );
  FloatType extent = l_high - l_low + r_high - r_low;
  return distance < extent;
}

template< typename Box >
bool
BetweenBoxes::operator ()(const Box& left, const Box& right) const
{
  typedef typename Box::value_type value_type;
  typedef typename Box::vector_type vector_type;

  const vector_type& l_low = left.low();
  const vector_type& l_high = left.high();
  const vector_type& r_low = right.low();
  const vector_type& r_high = right.high();

  return (
    one_dimensional_overlap< value_type >(
      l_low.elems[0],
      l_high.elems[0],
      r_low.elems[0],
      r_high.elems[0]
      )
    && one_dimensional_overlap< value_type >(
      l_low.elems[1],
      l_high.elems[1],
      r_low.elems[1],
      r_high.elems[1]
      )
    && one_dimensional_overlap< value_type >(
      l_low.elems[2],
      l_high.elems[2],
      r_low.elems[2],
      r_high.elems[2]
      )
    );
}

