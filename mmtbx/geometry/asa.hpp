#ifndef MMTBX_GEOMETRY_ASA_H
#define MMTBX_GEOMETRY_ASA_H

#include <mmtbx/geometry/primitive.hpp>

namespace mmtbx
{

namespace geometry
{

namespace asa
{

template< typename Vector >
class Sphere : public primitive::Sphere< Vector >
{
public:
  typedef typename primitive::Traits< Vector >::vector_type vector_type;
  typedef typename primitive::Traits< Vector >::value_type value_type;
  typedef size_t index_type;

private:
  size_t index_;

public:
  Sphere(const vector_type& centre, const value_type& radius, const size_t& index);
  ~Sphere();

  const index_type& index() const;
  vector_type low() const;
  vector_type high() const;
};

template< typename Vector >
class Transform
{
public:
  typedef Vector vector_type;
  typedef typename Vector::value_type value_type;
  typedef vector_type result_type;

private:
  vector_type centre_;
  value_type radius_;

public:
  Transform(const vector_type& centre, const value_type& radius);
  ~Transform();

  inline vector_type operator ()(const vector_type& point) const;
};

template< typename Object, typename Algorithm >
class OverlapEqualityFilter : private Algorithm
{
public:
  typedef Object object_type;
  typedef Algorithm algorithm_type;

private:
  object_type object_;

public:
  OverlapEqualityFilter(const object_type& object);
  ~OverlapEqualityFilter();

  inline bool operator ()(const object_type& other) const;
};

#include "asa.hxx"

} // namespace asa
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_ASA_H
