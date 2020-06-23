template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::Edge()
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::~Edge()
{}


template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::begin()
{
  node_type& n = this->node();
  return n.begin();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::end()
{
  node_type& n = this->node();
  return n.end();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::const_iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::begin() const
{
  node_type const& n = this->node();
  return n.begin();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::const_iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::end() const
{
  node_type const& n = this->node();
  return n.end();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::find(glyph_type const& key)
{
  node_type& n = this->node();
  return n.find( key );
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::const_iterator
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::find(glyph_type const& key) const
{
  node_type const& n = this->node();
  return n.find( key );
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
std::pair<
  typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::iterator,
  bool
  >
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::insert(value_type const& value)
{
  node_type& n = this->node();
  return n.insert( value );
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_parent() const
{
  return promote_link( this->parent() );
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_suffix() const
{
  return promote_link( this->suffix() );
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::promote_link(
  weak_ptr_type const& link
  ) const
{
  ptr_type locked = link.lock();

  if ( ! locked )
  {
    throw unavailable();
  }

  return locked;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
void
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::set_parent(
  ptr_type const& parent
  )
{
  this->parent() = parent;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_child_with_label(
  glyph_type const& label
  ) const
{
  const_iterator it = find( label );

  if ( it == end() )
  {
    throw nonexistent();
  }

  return it->second;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
void
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::attach_child(
  ptr_type const& child,
  glyph_type const& label
  )
{
  std::pair< typename edge_type::iterator, bool > res = insert( value_type( label, child ) );

  if ( ! res.second )
  {
    res.first->second = child;
  }
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::attach_child_if_not_present(
  ptr_type const& child,
  glyph_type const& label
  )
{
  return insert( value_type( label, child ) ).second;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::root()
{
  typedef Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > root_type;
  return boost::make_shared< root_type >();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::branch(
  index_type const& start,
  index_type const& stop
  )
{
  typedef Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > branch_type;
  return boost::make_shared< branch_type >(start, stop);
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::ptr_type
Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::leaf(
  index_type const& start,
  word_length_type const& word_length,
  suffix_label_type const& suffix_label
  )
{
  typedef Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > leaf_type;
  return boost::make_shared< leaf_type >(start, word_length, suffix_label);
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::Root()
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::~Root()
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::empty() const
{
  node_type const& n = this->node();
  return n.empty();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type 
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_start() const
{
  return index_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
void
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::set_start(
  index_type const& value
  )
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_stop() const
{
  return index_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix_label_type const&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::label() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type const&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node() const
{
  return node_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node()
{
  return node_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent()
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix()
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_root() const
{
  return true;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Root< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_leaf() const
{
  return false;
}


template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::Branch(
  index_type const& start,
  index_type const& stop
  )
  : start_( start ), stop_( stop )
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::~Branch()
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::empty() const
{
  node_type const& n = this->node();
  return n.empty();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_start() const
{
  return start_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
void
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::set_start(
  index_type const& value
  )
{
  start_ = value;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_stop() const
{
  return stop_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix_label_type const&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::label() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type const&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node() const
{
  return node_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node()
{
  return node_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent() const
{
  return parent_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent()
{
  return parent_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix() const
{
  return suffix_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix()
{
  return suffix_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_root() const
{
  return false;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Branch< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_leaf() const
{
  return false;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::Leaf(
  index_type const& start,
  word_length_type const& word_length,
  suffix_label_type const& suffix_label
  )
  : start_( start ), word_length_( word_length ), suffix_label_( suffix_label )
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::~Leaf()
{}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::empty() const
{
  return true;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_start() const
{
  return start_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
void
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::set_start(
  index_type const& value
  )
{
  start_ = value;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::index_type
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::get_stop() const
{
  return *word_length_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix_label_type const&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::label() const
{
  return suffix_label_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type const&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node_type&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::node()
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent() const
{
  return parent_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::parent()
{
  return parent_ptr_;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type const&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix() const
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
typename Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::weak_ptr_type&
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::suffix()
{
  throw bad_edge_type();
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_root() const
{
  return false;
}

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
bool
Leaf< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >::is_leaf() const
{
  return true;
}

