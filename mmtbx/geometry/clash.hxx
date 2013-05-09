template< typename Altloc >
bool
RegularAltlocStrategy< Altloc >::is_interacting_with(const strategy_type& other) const
{
  return true;
}

template< typename Altloc >
bool
RegularAltlocStrategy< Altloc >::is_interacting_with_alternate(Altloc identifier) const
{
  return true;
}

template< typename Altloc >
AlternateAltlocStrategy< Altloc >::AlternateAltlocStrategy(Altloc identifier)
  : identifier_( identifier )
{}

template< typename Altloc >
AlternateAltlocStrategy< Altloc >::~AlternateAltlocStrategy()
{}

template< typename Altloc >
bool
AlternateAltlocStrategy< Altloc >::is_interacting_with(const strategy_type& other) const
{
  return other.is_interacting_with_alternate( identifier_ );
}

template< typename Altloc >
bool
AlternateAltlocStrategy< Altloc >::is_interacting_with_alternate(
  Altloc identifier
  )
  const
{
  return identifier_ == identifier;
}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
Sphere< Vector, Identifier, Altloc, SymOp >::Sphere(
  const vector_type& centre,
  const value_type& radius,
  const identifier_type& molecule,
  const identifier_type& atom,
  const altloc_strategy_type& altloc_strategy,
  const symop_type& symop
  )
  : primitive::Sphere< Vector >( centre, radius ),
    molecule_( molecule ),
    atom_( atom ),
    altloc_strategy_( altloc_strategy ),
    symop_( symop )
{}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
Sphere< Vector, Identifier, Altloc, SymOp >::~Sphere()
{}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
const typename Sphere< Vector, Identifier, Altloc, SymOp >::identifier_type&
Sphere< Vector, Identifier, Altloc, SymOp >::molecule() const
{
  return molecule_;
}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
const typename Sphere< Vector, Identifier, Altloc, SymOp >::identifier_type&
Sphere< Vector, Identifier, Altloc, SymOp >::atom() const
{
  return atom_;
}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
const typename Sphere< Vector, Identifier, Altloc, SymOp >::altloc_strategy_type&
Sphere< Vector, Identifier, Altloc, SymOp >::altloc_strategy() const
{
  return altloc_strategy_;
}

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
const typename Sphere< Vector, Identifier, Altloc, SymOp >::symop_type&
Sphere< Vector, Identifier, Altloc, SymOp >::symop() const
{
  return symop_;
}

// Overlap and equivalence
template< typename Object, typename Algorithm >
OverlapInteractionFilter< Object, Algorithm >::OverlapInteractionFilter(
  const object_type& object,
  const value_type& tolerance
  )
  : object_( object ), tolerance_( tolerance )
{}

template< typename Object, typename Algorithm >
OverlapInteractionFilter< Object, Algorithm >::~OverlapInteractionFilter()
{}

template< typename Object, typename Algorithm >
bool
OverlapInteractionFilter< Object, Algorithm >::operator ()(
  const object_type& other
  )
  const
{
  return (
    Algorithm::operator ()( object_, other, tolerance_ )
    && (
      ( other.molecule() != object_.molecule() )
      || ( other.symop().t() != object_.symop().t() )
      || ( other.symop().r() != object_.symop().r() )
      )
    && object_.altloc_strategy()->is_interacting_with( *( other.altloc_strategy() ) )
    );
}
