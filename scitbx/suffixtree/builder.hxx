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
SuffixLinker< EdgePtr >::SuffixLinker()
{}

template< typename EdgePtr >
SuffixLinker< EdgePtr >::~SuffixLinker()
{}

template< typename EdgePtr >
void
SuffixLinker< EdgePtr >::operator ()(edge_ptr_type const& next)
{
  if ( previous_ )
  {
    previous_->suffix() = next;
  }

  previous_ = next;
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
  if ( not is_valid() )
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
  if ( not is_attached() )
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
      break;
    }
    catch ( nonexistent& e )
    {}

    catch ( mismatch& e )
    {
      position_.break_edge_here();
      suffix_linker( position_.get_edge_ptr() );
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
  suffix_linker( position_.get_edge_ptr() );
}

