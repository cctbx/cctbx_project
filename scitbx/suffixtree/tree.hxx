template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
Tree< Word, SuffixLabel, NodeAdaptor >::Tree()
  : root_( edge_type::root() ),
    word_ptr_( boost::make_shared< word_type >() ),
    construction_ptr_( boost::make_shared< bool >( false ) )
{}

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
Tree< Word, SuffixLabel, NodeAdaptor >::~Tree()
{}

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
typename Tree< Word, SuffixLabel, NodeAdaptor >::const_edge_ptr_type
Tree< Word, SuffixLabel, NodeAdaptor >::root() const
{
  return root_;
}

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
typename Tree< Word, SuffixLabel, NodeAdaptor >::word_type const&
Tree< Word, SuffixLabel, NodeAdaptor >::word() const
{
  return *word_ptr_;
}

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
bool
Tree< Word, SuffixLabel, NodeAdaptor >::in_construction() const
{
  return *construction_ptr_;
}


template< typename Tree >
typename Linkage< Tree >::edge_ptr_type
Linkage< Tree >::get_parent(edge_ptr_type const& edge_ptr)
{
  edge_ptr_type parent = ( edge_ptr->parent() ).lock();

  if ( not parent )
  {
    throw bad_state();
  }

  return parent;
}

template< typename Tree >
typename Linkage< Tree >::edge_ptr_type
Linkage< Tree >::get_suffix(edge_ptr_type const& edge_ptr)
{
  edge_ptr_type suffix = ( edge_ptr->suffix() ).lock();

  if ( not suffix )
  {
    throw bad_state();
  }

  return suffix;
}

template< typename Tree >
typename Linkage< Tree >::edge_ptr_type
Linkage< Tree >::get_child_with_label(
  edge_ptr_type const& edge_ptr,
  glyph_type const& glyph
  )
{
  typename edge_type::iterator it = edge_ptr->find( glyph );

  if ( it == edge_ptr->end() )
  {
    throw bad_state();
  }

  return it->second;
}


template< typename Tree >
typename Movement< Tree >::edge_ptr_type const&
Movement< Tree >::get_edge_ptr(cursor_type const& cursor)
{
  return cursor.first;
}

template< typename Tree >
typename Movement< Tree >::edge_ptr_type&
Movement< Tree >::get_edge_ptr(cursor_type& cursor)
{
  return cursor.first;
}

template< typename Tree >
typename Movement< Tree >::index_type&
Movement< Tree >::get_index(cursor_type& cursor)
{
  return cursor.second;
}

template< typename Tree >
typename Movement< Tree >::index_type const&
Movement< Tree >::get_index(cursor_type const& cursor)
{
  return cursor.second;
}

template< typename Tree >
bool
Movement< Tree >::is_at_edge_bottom(cursor_type const& cursor)
{
  return get_edge_ptr( cursor )->stop() == get_index( cursor );
}

template< typename Tree >
void
Movement< Tree >::set_to_edge_top(cursor_type& cursor, edge_ptr_type const& edge_ptr)
{
  get_edge_ptr( cursor ) = edge_ptr;
  get_index( cursor ) = edge_ptr->start();
}

template< typename Tree >
void
Movement< Tree >::forth(cursor_type& cursor)
{
  ++get_index( cursor );
}

template< typename Tree >
typename Movement< Tree >::cursor_type
Movement< Tree >::get_path_jump_destination(
  edge_ptr_type edge_ptr,
  word_iterator path_begin,
  word_iterator path_end
  )
{
  while ( true )
  {
    length_type edge_length = edge_ptr->stop() - edge_ptr->start();
    typename word_iterator::difference_type path_length = path_end - path_begin;

    if ( path_length <= edge_length )
    {
      return cursor_type( edge_ptr, edge_ptr->start() + path_length );
    }

    path_begin += edge_length;
    edge_ptr = linkage::get_child_with_label( edge_ptr, *path_begin );
  }
}

template< typename Tree >
typename Movement< Tree >::cursor_type
Movement< Tree >::get_suffix_position(cursor_type const& cursor, word_type const& word)
{
  edge_ptr_type const& current = get_edge_ptr( cursor );

  if ( not current->is_root() )
  {
    word_iterator path_begin = word.get_iterator_to( current->start() );

    edge_ptr_type suffix_ptr;
    edge_ptr_type parent_ptr = linkage::get_parent( current );

    if ( parent_ptr->is_root() )
    {
      suffix_ptr = parent_ptr;
      ++path_begin;
    }
    else
    {
      suffix_ptr = linkage::get_suffix( current );
    }

    word_iterator path_end = word.get_iterator_to( get_index( cursor ) );

    if ( path_begin == path_end )
    {
      return cursor_type( suffix_ptr, suffix_ptr->stop() );
    }
    else
    {
      return get_path_jump_destination(
        linkage::get_child_with_label( suffix_ptr, *path_begin ),
        path_begin,
        word.get_iterator_to( get_index( cursor ) )
        );
    }
  }
  else
  {
    return cursor;
  }
}

