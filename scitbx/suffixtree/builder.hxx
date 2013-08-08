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
    position_( tree_root_, const_cast< edge_type const& >( *tree_root_ ).start() ),
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
    if ( movement::is_at_edge_bottom( position_ ) )
    {
      edge_ptr_type& edge_ptr = movement::get_edge_ptr( position_ );
      typename edge_type::iterator it = edge_ptr->find( glyph );

      if ( it != edge_ptr->end() )
      {
        movement::set_to_edge_top( position_, it->second );
        movement::forth( position_ );
        break;
      }
    }

    else
    {

      if ( glyph == word[ movement::get_index( position_ ) ] )
      {
        movement::forth( position_ );
        break;
      }

      else
      {
        edge_ptr_type const& edge_ptr = movement::get_edge_ptr( position_ );
        index_type const& index = movement::get_index( position_ );
        index_type start = edge_ptr->start();
        edge_ptr_type new_branch_ptr = edge_type::branch( start, index );

        edge_weak_ptr_type parent_weak_ptr = edge_ptr->parent();
        new_branch_ptr->parent() = parent_weak_ptr;
        edge_ptr_type parent_ptr = parent_weak_ptr.lock();
        assert ( parent_ptr ); // root held
        typename edge_type::iterator pit = parent_ptr->find( word[ start ] );
        assert ( pit != parent_ptr->end() );
        pit->second = new_branch_ptr;

        edge_ptr->start() = position_.second;
        edge_ptr->parent() = new_branch_ptr;

        std::pair< typename edge_type::iterator, bool > res = new_branch_ptr->insert(
          typename edge_type::value_type( word[ index ], edge_ptr )
          );
        assert ( res.second );

        suffix_linker( new_branch_ptr );

        movement::get_edge_ptr( position_ ) = new_branch_ptr;
      }
    }
    edge_ptr_type& current = movement::get_edge_ptr( position_ );
    edge_ptr_type new_leaf = edge_type::leaf( phase_, word.length_ptr(), extension_ );
    std::pair< typename edge_type::iterator, bool > res = current->insert(
      typename edge_type::value_type( word[ phase_ ], new_leaf )
      );
    assert ( res.second );
    new_leaf->parent() = current;

    position_ = movement::get_suffix_position( position_, word );
    ++extension_;
  }

  ++phase_;
  suffix_linker( position_.first );
}

