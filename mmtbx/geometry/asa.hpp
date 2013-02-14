#ifndef MMTBX_GEOMETRY_ASA_H
#define MMTBX_GEOMETRY_ASA_H

#include <mmtbx/geometry/primitive.hpp>

#include <vector>

namespace mmtbx
{

namespace geometry
{

namespace asa
{

template< typename Vector >
class Sphere : public primitive::Sphere< Vector >
{
private:
  static size_t current;

public:
  typedef typename primitive::Traits< Vector >::vector_type vector_type;
  typedef typename primitive::Traits< Vector >::value_type value_type;
  typedef size_t index_type;

private:
  size_t index_;

public:
  Sphere(const vector_type& centre, const value_type& radius);
  ~Sphere();

  const index_type& index() const;
  vector_type low() const;
  vector_type high() const;
};

template< typename Vector >
bool
operator ==(const Sphere< Vector >& left, const Sphere< Vector >& right);

template< typename Vector >
bool
operator !=(const Sphere< Vector >& left, const Sphere< Vector >& right);

template< typename Sphere >
size_t
hash_value(const Sphere& object);


#include "asa.hxx"

} // namespace asa
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_ASA_H
