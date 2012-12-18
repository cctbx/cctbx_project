#ifndef MMTBX_GEOMETRY_PRIMITIVE_H
#define MMTBX_GEOMETRY_PRIMITIVE_H

namespace mmtbx
{

namespace geometry
{

namespace primitive
{

template< typename Vector >
struct Traits
{
public:
  typedef Vector vector_type;
  typedef typename Vector::value_type value_type;
};

template< typename Vector >
class Sphere
{
public:
  typedef typename Traits< Vector >::vector_type vector_type;
  typedef typename Traits< Vector >::value_type value_type;

private:
  vector_type centre_;
  value_type radius_;
  value_type radius_sq_;

public:
  Sphere(const vector_type& centre, const value_type& radius);
  ~Sphere();

  inline const vector_type& centre() const;
  inline const value_type& radius() const;
  inline const value_type& radius_sq() const;
};

template< typename Vector >
class Box
{
public:
  typedef typename Traits< Vector >::vector_type vector_type;
  typedef typename Traits< Vector >::value_type value_type;

private:
  vector_type low_;
  vector_type high_;

  // Private for quick construction
  Box(const vector_type& low, const vector_type& high);

public:
  ~Box();

  inline const vector_type& low() const;
  inline const vector_type& high() const;

  static Box< Vector > from_corners(
    const vector_type& corner1,
    const vector_type& corner2
    );
  static Box< Vector > around_sphere(
    const vector_type& centre,
    const value_type& radius
    );
  static Box< Vector > from_sphere(const Sphere< Vector >& sphere);
};

template< typename Vector >
class BSphere : public Sphere< Vector >, public Box< Vector >
{
public:
  typedef typename Traits< Vector >::vector_type vector_type;
  typedef typename Traits< Vector >::value_type value_type;

public:
  BSphere(const vector_type& centre, const value_type& radius);
  ~BSphere();
};


#include "primitive.hxx"

} // namespace primitive
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_PRIMITIVE_H
