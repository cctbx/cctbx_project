template< typename Value >
Character< Value >::Character(const value_type& value)
  : value_( value )
{}

template< typename Value >
Character< Value >::~Character()
{}

template< typename Value >
const typename Character< Value >::value_type&
Character< Value >::value() const
{
  return value_;
}

template< typename Value >
std::ostream&
operator <<(std::ostream& stream, const Character< Value >& character)
{
  return stream << character.value();
}

template< typename Value >
std::size_t
hash_value(const Character< Value >& character)
{
  return hash_value( character.value() );
}

template< typename Value >
bool
operator ==(const Character< Value >& lhs, const Character< Value >& rhs)
{
  return lhs.value() == rhs.value();
}

template< typename Value >
bool
operator <(const Character< Value >& lhs, const Character< Value >& rhs)
{
  return lhs.value() < rhs.value();
}


template< typename Index >
Terminator< Index >::Terminator(const index_type& index)
  : index_( index )
{}

template< typename Index >
Terminator< Index >::~Terminator()
{}

template< typename Index >
const typename Terminator< Index >::index_type&
Terminator< Index >::index() const
{
  return index_;
}

template< typename Index >
std::ostream& operator <<(std::ostream& stream, const Terminator< Index >& terminator)
{
  return stream << "T(" << terminator.index() << ")";
}

template< typename Index >
std::size_t
hash_value(const Terminator< Index >& terminator)
{
  return hash_value( terminator.index() );
}

template< typename Index >
bool
operator ==(const Terminator< Index >& lhs, const Terminator< Index >& rhs)
{
  return lhs.index() == rhs.index();
}

template< typename Index >
bool
operator <(const Terminator< Index >& lhs, const Terminator< Index >& rhs)
{
  return lhs.index() < rhs.index();
}
