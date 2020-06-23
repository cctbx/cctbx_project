#ifndef SUFFIXTREE_MATCHING_STATISTICS_HPP_
#define SUFFIXTREE_MATCHING_STATISTICS_HPP_

#include <scitbx/suffixtree/tree.hpp>
#include <scitbx/suffixtree/exception.hpp>

#include <iterator>
#include <utility>

namespace scitbx
{

namespace suffixtree
{

template< typename Tree, typename InputIterator >
class MSI : public std::iterator<
  std::input_iterator_tag,
  std::pair<
    typename Tree::length_type,
    std::pair< typename Tree::edge_ptr_type, typename Tree::index_type >
    > const
  >
{
public:
  typedef Tree tree_type;
  typedef typename tree_type::length_type length_type;
  typedef typename tree_type::edge_ptr_type edge_ptr_type;
  typedef typename tree_type::index_type index_type;
  typedef typename tree_type::word_type word_type;
  typedef typename tree_type::glyph_type glyph_type;
  typedef typename tree_type::cursor_type cursor_type;

  typedef std::pair< edge_ptr_type, index_type > position_type;
  typedef std::pair< length_type, position_type > value_type;
  typedef value_type const* pointer_type;
  typedef value_type const& reference_type;

  typedef InputIterator text_iterator_type;

  typedef MSI iterator;

private:
  cursor_type cursor_;

  text_iterator_type begin_;
  text_iterator_type end_;
  length_type matching_;

  value_type result_;

public:
  MSI(
    tree_type const& tree,
    text_iterator_type begin,
    text_iterator_type end
    );
  ~MSI();

  reference_type operator *() const;
  pointer_type operator ->() const;

  iterator& operator ++();
  iterator operator ++(int);

private:
  void follow_until_mismatch();

  template< typename T, typename II >
  friend bool operator ==(MSI< T, II > const& lhs, MSI< T, II > const& rhs);
};

template< typename T, typename II >
bool operator !=(MSI< T, II > const& lhs, MSI< T, II > const& rhs);


// Definitions
template< typename Tree, typename InputIterator >
MSI< Tree, InputIterator >::MSI(
  tree_type const& tree,
  text_iterator_type begin,
  text_iterator_type end
  )
 : cursor_( tree.cursor() ),
   begin_( begin ), end_( end ), matching_( 0 )
{
  if ( tree.in_construction() )
  {
    throw bad_tree();
  }

  follow_until_mismatch();
}

template< typename Tree, typename InputIterator >
MSI< Tree, InputIterator >::~MSI()
{}

template< typename Tree, typename InputIterator >
typename MSI< Tree, InputIterator >::iterator&
MSI< Tree, InputIterator >::operator ++()
{
  cursor_.to_suffix_position();

  if ( 0 < matching_ )
  {
    --matching_;
  }
  else
  {
    ++begin_;
  }

  follow_until_mismatch();

  return *this;
}

template< typename Tree, typename InputIterator >
typename MSI< Tree, InputIterator >::iterator
MSI< Tree, InputIterator >::operator ++(int)
{
  iterator old( *this );
  ++( *this );
  return old;
}


template< typename Tree, typename InputIterator >
typename MSI< Tree, InputIterator >::reference_type
MSI< Tree, InputIterator >::operator *() const
{
  return result_;
}

template< typename Tree, typename InputIterator >
typename MSI< Tree, InputIterator >::pointer_type
MSI< Tree, InputIterator >::operator ->() const
{
  return &result_;
}

template< typename Tree, typename InputIterator >
void
MSI< Tree, InputIterator >::follow_until_mismatch()
{
  while( begin_ != end_ )
  {
    try
    {
      cursor_.forth_with( *begin_ );
    }
    catch ( nonexistent& e )
    {
      break;
    }
    catch ( mismatch& e )
    {
      break;
    }

    ++matching_;
    ++begin_;
  }

  result_.first = matching_;
  result_.second.first = cursor_.get_edge_ptr();
  result_.second.second = cursor_.get_index();
}

template< typename Tree, typename InputIterator >
bool operator ==(
  MSI< Tree, InputIterator > const& lhs,
  MSI< Tree, InputIterator > const& rhs
  )
{
  return ( lhs.begin_ == rhs.begin_
    && lhs.end_ == rhs.end_
    && lhs.matching_ == rhs.matching_ );
}

template< typename Tree, typename InputIterator >
bool operator !=(
  MSI< Tree, InputIterator > const& lhs,
  MSI< Tree, InputIterator > const& rhs
  )
{
  return ! ( lhs == rhs );
}




} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_MATCHING_STATISTICS_HPP_
