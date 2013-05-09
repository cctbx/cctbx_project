#ifndef MMTBX_GEOMETRY_CLASH_H
#define MMTBX_GEOMETRY_CLASH_H

#include <mmtbx/geometry/primitive.hpp>

#include <boost/shared_ptr.hpp>

namespace mmtbx
{

namespace geometry
{

namespace clash
{

template< typename Altloc >
class AltlocStrategy
{
public:
  typedef AltlocStrategy< Altloc > strategy_type;

public:
  virtual bool is_interacting_with(const strategy_type& other) const = 0;
  virtual bool is_interacting_with_alternate(Altloc identifier) const = 0;
};

template< typename Altloc >
class RegularAltlocStrategy : public AltlocStrategy< Altloc >
{
public:
  typedef AltlocStrategy< Altloc > strategy_type;

public:
  virtual bool is_interacting_with(const strategy_type& other) const;
  virtual bool is_interacting_with_alternate(Altloc identifier) const;
};

template< typename Altloc >
class AlternateAltlocStrategy : public AltlocStrategy< Altloc >
{
public:
  typedef AltlocStrategy< Altloc > strategy_type;

private:
  Altloc identifier_;

public:
  AlternateAltlocStrategy(Altloc identifier);
  ~AlternateAltlocStrategy();

  virtual bool is_interacting_with(const strategy_type& other) const;
  virtual bool is_interacting_with_alternate(Altloc identifier) const;
};

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
class Sphere : public primitive::Sphere< Vector >
{
public:
  typedef typename primitive::Traits< Vector >::vector_type vector_type;
  typedef typename primitive::Traits< Vector >::value_type value_type;
  typedef Identifier identifier_type;
  typedef Altloc altloc_type;
  typedef boost::shared_ptr< AltlocStrategy< altloc_type > > altloc_strategy_type;
  typedef SymOp symop_type;

private:
   identifier_type molecule_;
   identifier_type atom_;
   altloc_strategy_type altloc_strategy_;
   symop_type symop_;

public:
  Sphere(
    const vector_type& centre,
    const value_type& radius,
    const identifier_type& molecule,
    const identifier_type& atom,
    const altloc_strategy_type& altloc_strategy,
    const symop_type& symop
    );
  ~Sphere();

  const identifier_type& molecule() const;
  const identifier_type& atom() const;
  const altloc_strategy_type& altloc_strategy() const;
  const symop_type& symop() const;
};

template< typename Object, typename Algorithm >
class OverlapInteractionFilter : private Algorithm
{
public:
  typedef Object object_type;
  typedef typename Object::value_type value_type;

private:
  object_type object_;
  value_type tolerance_;

public:
  OverlapInteractionFilter(const object_type& object, const value_type& tolerance);
  ~OverlapInteractionFilter();

  inline bool operator ()(const object_type& other) const;
};

#include "clash.hxx"

} // namespace clash
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_CLASH_H
