#ifndef SUFFIXTREE_GLYPH_HPP_
#define SUFFIXTREE_GLYPH_HPP_

#include <boost/functional/hash.hpp>
#include <ostream>

namespace scitbx
{

namespace suffixtree
{

namespace glyph
{

template< typename Value >
class Character
{
public:
  typedef Value value_type;

private:
  value_type value_;

public:
  explicit Character(const value_type& value);
  ~Character();

  const value_type& value() const;
};

template< typename Value >
std::ostream& operator <<(std::ostream& stream, const Character< Value >& character);

template< typename Value >
std::size_t hash_value(const Character< Value >& character);

template< typename Value >
bool operator ==(const Character< Value >& lhs, const Character< Value >& rhs);

template< typename Value >
bool operator <(const Character< Value >& lhs, const Character< Value >& rhs);

template< typename Index >
class Terminator
{
public:
  typedef Index index_type;

private:
  index_type index_;

public:
  explicit Terminator(const index_type& index);
  ~Terminator();

  const index_type& index() const;
};

template< typename Index >
std::ostream& operator <<(std::ostream& stream, const Terminator< Index >& terminator);

template< typename Index >
std::size_t hash_value(const Terminator< Index >& character);

template< typename Index >
bool operator <(const Terminator< Index >& lhs, const Terminator< Index >& rhs);

template< typename Index >
bool operator ==(const Terminator< Index >& lhs, const Terminator< Index >& rhs);

#include "glyph.hxx"

} // namespace glyph
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_GLYPH_HPP_
