template< typename Tree >
Builder< Tree >::Builder()
{}

template< typename Tree >
Builder< Tree >::~Builder()
{}

template< typename Tree >
typename Builder< Tree >::edge_ptr_type const&
Builder< Tree >::tree_root(tree_type const& tree)
{
  return tree.root_;
}

template< typename Tree >
typename Builder< Tree >::word_ptr_type const&
Builder< Tree >::tree_word_ptr(tree_type const& tree)
{
  return tree.word_ptr_;
}

template< typename Tree >
typename Builder< Tree >::construction_ptr_type const&
Builder< Tree >::tree_construction_ptr(tree_type const& tree)
{
  return tree.construction_ptr_;
}

template< typename EdgePtr >
SuffixLinkerState< EdgePtr >::~SuffixLinkerState()
{}

template< typename EdgePtr >
bool
SuffixLinkerEmpty< EdgePtr >::process_existing(
  edge_ptr_type& stored,
  edge_ptr_type const& next
  ) const
{
  return false;
}

template< typename EdgePtr >
bool
SuffixLinkerEmpty< EdgePtr >::process_new(
  edge_ptr_type& stored,
  edge_ptr_type const& next
  ) const
{
  stored = next;
  return true;
}

template< typename EdgePtr >
bool
SuffixLinkerPrimed< EdgePtr >::process_existing(
  edge_ptr_type& stored,
  edge_ptr_type const& next
  ) const
{
  stored->suffix() = next;
  return true;
}

template< typename EdgePtr >
bool
SuffixLinkerPrimed< EdgePtr >::process_new(
  edge_ptr_type& stored,
  edge_ptr_type const& next
  ) const
{
  stored->suffix() = next;
  stored = next;
  return false;
}

template< typename EdgePtr >
typename SuffixLinker< EdgePtr >::empty_state_type const
SuffixLinker< EdgePtr >::empty_state = empty_state_type();

template< typename EdgePtr >
typename SuffixLinker< EdgePtr >::primed_state_type const
SuffixLinker< EdgePtr >::primed_state = primed_state_type();

template< typename EdgePtr >
SuffixLinker< EdgePtr >::SuffixLinker()
  : state_ ( &empty_state )
{}

template< typename EdgePtr >
SuffixLinker< EdgePtr >::~SuffixLinker()
{}

template< typename EdgePtr >
void
SuffixLinker< EdgePtr >::process_existing(edge_ptr_type const& next)
{
  if ( state_->process_existing( storage_, next ) )
  {
    state_ = &empty_state;
  }
}

template< typename EdgePtr >
void
SuffixLinker< EdgePtr >::process_new(edge_ptr_type const& next)
{
  if ( state_->process_new( storage_, next ) )
  {
    state_ = &primed_state;
  }
}

template< typename Tree >
Ukkonen< Tree >::Ukkonen(tree_type const& tree)
  : tree_root_( Builder< Tree >::tree_root( tree ) ),
    tree_word_ptr_( Builder< Tree >::tree_word_ptr( tree ) ),
    tree_construction_ptr_( Builder< Tree >::tree_construction_ptr( tree ) ),
    position_( tree_root_, tree_word_ptr_ ),
    phase_( *( tree_word_ptr_->length_ptr() ) ),
    extension_( phase_ ),
    is_attached_( true )
{
  if ( *tree_construction_ptr_ )
  {
    throw bad_state();
  }

  *tree_construction_ptr_ = true;
}

template< typename Tree >
Ukkonen< Tree >::~Ukkonen()
{}

template< typename Tree >
bool
Ukkonen< Tree >::is_attached() const
{
  return is_attached_;
}

template< typename Tree >
bool
Ukkonen< Tree >::is_valid() const
{
  return phase_ == extension_;
}

template< typename Tree >
void
Ukkonen< Tree >::detach()
{
  if ( ! is_valid() )
  {
    throw bad_state();
  }

  tree_root_.reset();
  tree_word_ptr_.reset();

  *tree_construction_ptr_ = false;
  tree_construction_ptr_.reset();

  is_attached_ = false;
}


template< typename Tree >
void
Ukkonen< Tree >::push_back(glyph_type const& glyph)
{
  if ( ! is_attached() )
  {
    throw bad_state();
  }

  word_type& word = *tree_word_ptr_;
  word.push_back( glyph );
  suffix_linker_type suffix_linker;

  while ( extension_ <= phase_ )
  {
    try
    {
      position_.forth_with( glyph );
      suffix_linker.process_existing( position_.get_edge_ptr()->get_parent() );
      break;
    }
    catch ( nonexistent& e )
    {
      suffix_linker.process_existing( position_.get_edge_ptr() );
    }

    catch ( mismatch& e )
    {
      position_.break_edge_here();
      suffix_linker.process_new( position_.get_edge_ptr() );
    }

    edge_ptr_type const& current = position_.get_edge_ptr();
    edge_ptr_type new_leaf = edge_type::leaf( phase_, word.length_ptr(), extension_ );
    bool res = current->attach_child_if_not_present( new_leaf, word[ phase_ ] );
    assert ( res );
    new_leaf->parent() = current;
    position_.to_suffix_position();

    ++extension_;
  }

  ++phase_;
}

