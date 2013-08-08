template< typename Glyph >
Single< Glyph >::Single()
  : length_ptr_( new index_type() )
{}

template< typename Glyph >
Single< Glyph >::~Single()
{}

template< typename Glyph >
void
Single< Glyph >::push_back(const glyph_type& glyph)
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
  return boost::const_pointer_cast< const index_type >( length_ptr_ );
}

template< typename Glyph >
typename Single< Glyph >::substring_type
Single< Glyph >::substring(const index_type& begin, const index_type& end)
  const
{
  return substring_type(  get_iterator_to( begin ), get_iterator_to( end ) );
}

template< typename Glyph >
typename Single< Glyph >::const_iterator
Single< Glyph >::get_iterator_to(const index_type& index)
  const
{
  const_iterator cit = data_.begin();
  std::advance( cit, index );
  return cit;
}

template< typename Glyph >
const typename Single< Glyph >::glyph_type&
Single< Glyph >::operator [](const index_type& index) const
{
  return data_[ index ];
}

template< typename Glyph >
typename Single< Glyph >::glyph_type&
Single< Glyph >::operator [](const index_type& index)
{
  return data_[ index ];
}


template< typename Glyph, typename Traits >
Multiple< Glyph, Traits >::Multiple()
  : length_ptr_( new index_type() )
{}

template< typename Glyph, typename Traits >
Multiple< Glyph, Traits >::~Multiple()
{}

template< typename Glyph, typename Traits >
void
Multiple< Glyph, Traits >::push_back(const glyph_type& glyph)
{
  data_.push_back( glyph );
  ++( *length_ptr_ );

  if (Traits::is_terminator( glyph ) )
  {
    set_word_boundary();
  }
}

template< typename Glyph, typename Traits >
typename Multiple< Glyph, Traits >::length_type
Multiple< Glyph, Traits >::length() const
{
  return boost::const_pointer_cast< const index_type >( length_ptr_ );
}

template< typename Glyph, typename Traits >
typename Multiple< Glyph, Traits >::substring_type
Multiple< Glyph, Traits >::substring(const index_type& begin, const index_type& end)
  const
{
  const_iterator b = data_.begin();
  std::advance( b, begin );

  const_iterator e = data_.begin();
  std::advance( e, end );

  return substring_type( b, e );
}

template< typename Glyph, typename Traits >
const typename Multiple< Glyph, Traits >::glyph_type&
Multiple< Glyph, Traits >::operator [](const index_type& index) const
{
  return data_[ index ];
}

template< typename Glyph, typename Traits >
typename Multiple< Glyph, Traits >::glyph_type&
Multiple< Glyph, Traits >::operator [](const index_type& index)
{
  return data_[ index ];
}

template< typename Glyph, typename Traits >
void
Multiple< Glyph, Traits >::set_word_boundary()
{
  length_ptr_ = boost::make_shared< index_type >( *length_ptr_ );
}

