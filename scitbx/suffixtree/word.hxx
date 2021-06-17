template< typename Glyph >
Single< Glyph >::Single()
  : length_ptr_( new length_type() )
{}

template< typename Glyph >
Single< Glyph >::~Single()
{}

template< typename Glyph >
void
Single< Glyph >::push_back(glyph_type const& glyph)
{
  data_.push_back( glyph );
  ++( *length_ptr_ );
}

template< typename Glyph >
typename Single< Glyph >::length_type
Single< Glyph >::size() const
{
  return *length_ptr_;
}

template< typename Glyph >
typename Single< Glyph >::const_length_ptr_type
Single< Glyph >::length_ptr() const
{
  return boost::const_pointer_cast< const_length_type >( length_ptr_ );
}

template< typename Glyph >
typename Single< Glyph >::substring_type
Single< Glyph >::substring(index_type const& begin, index_type const& end)
  const
{
  return substring_type(  get_iterator_to( begin ), get_iterator_to( end ) );
}

template< typename Glyph >
typename Single< Glyph >::const_iterator
Single< Glyph >::get_iterator_to(index_type const& index)
  const
{
  const_iterator cit = data_.begin();
  std::advance( cit, index );
  return cit;
}

template< typename Glyph >
typename Single< Glyph >::glyph_type const&
Single< Glyph >::operator [](index_type const& index) const
{
  return data_[ index ];
}

template< typename Glyph >
typename Single< Glyph >::glyph_type&
Single< Glyph >::operator [](index_type const& index)
{
  return data_[ index ];
}

