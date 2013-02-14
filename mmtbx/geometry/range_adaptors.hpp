#ifndef MMTBX_GEOMETRY_RANGE_ADAPTORS_H
#define MMTBX_GEOMETRY_RANGE_ADAPTORS_H

namespace mmtbx
{

namespace geometry
{

namespace adaptors
{

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

#include "range_adaptors.hxx"

} // namespace adaptors
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_RANGE_ADAPTORS_H
