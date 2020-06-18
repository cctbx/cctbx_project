template< typename Edge, typename Word >
Cursor< Edge, Word >::Cursor(edge_ptr_type const& edge_ptr, word_ptr_type const& word_ptr)
  : word_ptr_( word_ptr ), edge_ptr_( edge_ptr ),
    index_( edge_ptr_->get_start() )
{}

template< typename Edge, typename Word >
Cursor< Edge, Word >::~Cursor()
{}


template< typename Edge, typename Word >
typename Cursor< Edge, Word >::edge_ptr_type const&
Cursor< Edge, Word >::get_edge_ptr() const
{
  return edge_ptr_;
}

template< typename Edge, typename Word >
typename Cursor< Edge, Word >::index_type const&
Cursor< Edge, Word >::get_index() const
{
  return index_;
}

template< typename Edge, typename Word >
bool
Cursor< Edge, Word >::is_at_edge_top() const
{
  return get_edge_ptr()->get_start() == get_index();
}

template< typename Edge, typename Word >
bool
Cursor< Edge, Word >::is_at_edge_bottom() const
{
  return get_edge_ptr()->get_stop() == get_index();
}

template< typename Edge, typename Word >
bool
Cursor< Edge, Word >::is_at_leaf_bottom() const
{
  edge_ptr_type const& edge_ptr = get_edge_ptr();

  return ( edge_ptr->is_leaf() && edge_ptr->get_stop() == get_index() );
}

template< typename Edge, typename Word >
typename Cursor< Edge, Word >::glyph_type const&
Cursor< Edge, Word >::get_current_character() const
{
  return ( *word_ptr_ )[ get_index() ];
}

template< typename Edge, typename Word >
typename Cursor< Edge, Word >::word_type const&
Cursor< Edge, Word >::get_word() const
{
  return *word_ptr_;
}

/*
template< typename Edge, typename Word >
typename Cursor< Edge, Word >::edge_ptr_type
Cursor< Edge, Word >::get_child_with_label(
  cursor_type const& cursor,
  glyph_type const& label
  )
{
  return get_edge_ptr( cursor )->get_child_with_label( label );
}
*/

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::forth_to_child(glyph_type const& label)
{
  edge_ptr_ = get_edge_ptr()->get_child_with_label( label );
  index_ = edge_ptr_->get_start();
}

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::forth_on_edge()
{
  if ( is_at_edge_bottom() )
  {
    throw bad_state();
  }

  ++index_;
}

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::forth_with(glyph_type const& label)
{
  if ( is_at_edge_bottom() )
  {
    forth_to_child( label );
    forth_on_edge();
  }
  else
  {
    if ( get_current_character() == label )
    {
      forth_on_edge();
    }
    else
    {
      throw mismatch();
    }
  }
}

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::break_edge_here()
{
  if ( is_at_edge_bottom() || is_at_edge_top() )
  {
    throw bad_state();
  }

  index_type start = edge_ptr_->get_start();
  edge_ptr_type new_branch_ptr = edge_type::branch( start, index_ );

  edge_weak_ptr_type parent_weak_ptr = edge_ptr_->parent();
  new_branch_ptr->parent() = parent_weak_ptr;
  edge_ptr_type parent_ptr = parent_weak_ptr.lock();

  if ( ! parent_ptr )
  {
    throw bad_tree();
  }

  word_type const& word = *word_ptr_;
  typename edge_type::iterator pit = parent_ptr->find( word[ start ] );
  assert ( pit != parent_ptr->end() );
  pit->second = new_branch_ptr;

  edge_ptr_->set_start( index_ );
  edge_ptr_->parent() = new_branch_ptr;

  bool res = new_branch_ptr->attach_child_if_not_present( edge_ptr_, word[ index_ ] );
  assert ( res );

  edge_ptr_ = new_branch_ptr;
}

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::to_suffix_position()
{
  word_type const& word = *word_ptr_;

  if ( ! edge_ptr_->is_root() )
  {
    edge_ptr_type suffix_ptr;
    edge_ptr_type parent_ptr = edge_ptr_->get_parent();
    word_iterator path_begin = word.get_iterator_to( edge_ptr_->get_start() );

    if ( parent_ptr->is_root() )
    {
      suffix_ptr = parent_ptr;
      ++path_begin;
    }
    else
    {
      suffix_ptr = parent_ptr->get_suffix();
    }

    word_iterator path_end = word.get_iterator_to( index_ );

    if ( path_begin == path_end )
    {
      edge_ptr_ = suffix_ptr;
      index_ = suffix_ptr->get_stop();
    }
    else
    {
      path_jump_from_top_of(
        suffix_ptr->get_child_with_label( *path_begin ),
        path_begin,
        path_end
        );
    }
  }
}

template< typename Edge, typename Word >
void
Cursor< Edge, Word >::path_jump_from_top_of(
  edge_ptr_type edge_ptr,
  word_iterator path_begin,
  word_iterator path_end
  )
{
  while ( true )
  {
    length_type edge_length = edge_ptr->get_stop() - edge_ptr->get_start();
    typename word_iterator::difference_type path_length = path_end - path_begin;

    if ( path_length <= edge_length )
    {
      edge_ptr_ = edge_ptr;
      index_ = edge_ptr->get_start() + path_length;
      break;
    }

    path_begin += edge_length;
    edge_ptr = edge_ptr->get_child_with_label( *path_begin );
  }
}


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
typename Tree< Word, SuffixLabel, NodeAdaptor >::edge_ptr_type
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

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
typename Tree< Word, SuffixLabel, NodeAdaptor >::cursor_type
Tree< Word, SuffixLabel, NodeAdaptor >::cursor() const
{
  return cursor_type( root_, word_ptr_ );
}


template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
bool
operator ==(
  Tree< Word, SuffixLabel, NodeAdaptor > const& lhs,
  Tree< Word, SuffixLabel, NodeAdaptor > const& rhs
  )
{
  return ( lhs.root_ == rhs.root_
    && lhs.word_ptr_ == rhs.word_ptr
    && lhs.construction_ptr_ == rhs.construction_ptr_ );
}

