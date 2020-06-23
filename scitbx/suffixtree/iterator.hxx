template< typename Edge >
PreOrder< Edge >::PreOrder(ptr_type const& root, bool at_top)
  : root_( root ), at_top_( at_top )
{
  if ( ! root_->empty() )
  {
    pos_ = at_top ? root->begin() : root->end();
  }
}

template< typename Edge >
PreOrder< Edge >::~PreOrder()
{}

template< typename Edge >
typename PreOrder< Edge >::iterator&
PreOrder< Edge >::operator ++()
{
  if ( at_top_ )
  {
    at_top_ = false;
    return *this;
  }

  if ( ! ( pos_->second )->empty() )
  {
    underlying_iterators_deque_.push_back( pos_ );
    pos_ = ( pos_->second )->begin();
  }
  else
  {
    ++pos_;

    while ( !underlying_iterators_deque_.empty()
      && pos_ == ( ( underlying_iterators_deque_.back() )->second )->end() )
    {
      pos_ = underlying_iterators_deque_.back();
      underlying_iterators_deque_.pop_back();
      ++pos_;
    }
  }

  return *this;
}

template< typename Edge >
typename PreOrder< Edge >::iterator
PreOrder< Edge >::operator ++(int)
{
  iterator old( *this );
  ++( *this );
  return old;
}


template< typename Edge >
typename PreOrder< Edge >::reference
PreOrder< Edge >::operator *() const
{
  return at_top_ ? root_ : pos_->second;
}

template< typename Edge >
typename PreOrder< Edge >::pointer
PreOrder< Edge >::operator ->() const
{
  return at_top_ ? &root_ : &( pos_->second );
}

template< typename Edge >
typename PreOrder< Edge >::iterator
PreOrder< Edge >::begin(ptr_type const& root)
{
  return iterator( root, true );
}

template< typename Edge >
typename PreOrder< Edge >::iterator
PreOrder< Edge >::end(ptr_type const& root)
{
  return iterator( root, false );
}

template< typename Edge >
bool operator ==(PreOrder<Edge> const& lhs, PreOrder<Edge> const& rhs)
{
  return ( ( lhs.pos_ == rhs.pos_ ) && ( lhs.root_ == rhs.root_ )
    && ( lhs.at_top_ == rhs.at_top_ ) );
}

template< typename Edge >
bool operator !=(PreOrder<Edge> const& lhs, PreOrder<Edge> const& rhs)
{
  return ! ( lhs == rhs );
}


template< typename Edge >
PostOrder< Edge >::PostOrder(ptr_type const& root, bool at_top)
  : root_( root )
{
  if ( root_->empty() )
  {
    at_top_ = at_top;
  }
  else if ( ! at_top )
  {
    at_top_ = false;
    pos_ = root_->end();
  }
  else
  {
    at_top_ = false;
    pos_ = root_->begin();
    descend();
  }
}

template< typename Edge >
PostOrder< Edge >::~PostOrder()
{}

template< typename Edge >
typename PostOrder< Edge >::iterator&
PostOrder< Edge >::operator ++()
{
  if ( at_top_ )
  {
    at_top_ = false;
    return *this;
  }

  ++pos_;

  if ( underlying_iterators_deque_.empty() )
  {
    if ( pos_ == root_->end() )
    {
      at_top_ = true;
    }
    else
    {
      descend();
    }
  }
  else
  {
    const underlying_const_iterator& parent = underlying_iterators_deque_.back();

    if ( pos_ != ( parent->second )->end() )
    {
      descend();
    }
    else
    {
      pos_ = parent;
      underlying_iterators_deque_.pop_back();
    }
  }

  return *this;
}

template< typename Edge >
typename PostOrder< Edge >::iterator
PostOrder< Edge >::operator ++(int)
{
  iterator old( *this );
  ++( *this );
  return old;
}


template< typename Edge >
typename PostOrder< Edge >::reference
PostOrder< Edge >::operator *() const
{
  return at_top_ ? root_ : pos_->second;
}

template< typename Edge >
typename PostOrder< Edge >::pointer
PostOrder< Edge >::operator ->() const
{
  return at_top_ ? &root_ : &( pos_->second );
}

template< typename Edge >
typename PostOrder< Edge >::iterator
PostOrder< Edge >::begin(ptr_type const& root)
{
  return iterator( root, true );
}

template< typename Edge >
typename PostOrder< Edge >::iterator
PostOrder< Edge >::end(ptr_type const& root)
{
  return iterator( root, false );
}

template< typename Edge >
void
PostOrder< Edge >::descend()
{
  while ( ! ( pos_->second )->empty() )
  {
    underlying_iterators_deque_.push_back( pos_ );
    pos_ = ( pos_->second )->begin();
  }
}

template< typename Edge >
bool operator ==(PostOrder<Edge> const& lhs, PostOrder<Edge> const& rhs)
{
  return ( ( lhs.pos_ == rhs.pos_ ) && ( lhs.root_ == rhs.root_ )
    && ( lhs.at_top_ == rhs.at_top_ ) );
}

template< typename Edge >
bool operator !=(PostOrder<Edge> const& lhs, PostOrder<Edge> const& rhs)
{
  return ! ( lhs == rhs );
}
