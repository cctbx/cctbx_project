#ifndef SUFFIXTREE_ITERATOR_HPP_
#define SUFFIXTREE_ITERATOR_HPP_

#include <scitbx/suffixtree/edge.hpp>

#include <boost/mpl/if.hpp>

#include <iterator>
#include <deque>

namespace scitbx
{

namespace suffixtree
{

namespace iterator
{

template< typename Edge >
struct Selector
{
  typedef edge::Traits< Edge > traits;
  typedef typename traits::iterator iterator;
  typedef typename traits::ptr_type value_type;
  typedef value_type ptr_type;

  typedef typename boost::mpl::if_<
    typename traits::selector_type,
    typename Edge::const_ptr_type const&,
    typename Edge::ptr_type&
    >::type reference;

  typedef typename boost::mpl::if_<
    typename traits::selector_type,
    typename Edge::const_ptr_type* const,
    typename Edge::ptr_type*
    >::type pointer;
};

template< typename Edge >
class PreOrder : public std::iterator<
  std::forward_iterator_tag,
  typename Selector< Edge >::value_type,
  std::ptrdiff_t,
  typename Selector< Edge >::reference,
  typename Selector< Edge >::pointer
  >
{
public:
  typedef Edge edge_type;
  typedef Selector< Edge > selector_type;
  typedef typename selector_type::value_type value_type;
  typedef typename selector_type::reference reference_type;
  typedef typename selector_type::pointer pointer_type;
  typedef typename selector_type::iterator underlying_iterator;
  typedef typename selector_type::ptr_type ptr_type;

  typedef PreOrder iterator;

private:
  ptr_type root_;
  bool at_top_;
  underlying_iterator pos_;
  std::deque< underlying_iterator > underlying_iterators_deque_;

private:
  PreOrder(ptr_type const& root, bool at_top);

public:
  static iterator begin(ptr_type const& root);
  static iterator end(ptr_type const& root);

public:
  ~PreOrder();

  reference_type operator *();
  pointer_type operator ->();

  iterator& operator ++();
  iterator operator ++(int);

  template< typename E2 >
  friend bool operator ==(PreOrder< E2 > const& lhs, PreOrder< E2 > const& rhs);
};

template< typename Edge >
bool operator !=(PreOrder< Edge > const& lhs, PreOrder< Edge > const& rhs);

template< typename Edge >
class PostOrder : public std::iterator<
  std::forward_iterator_tag,
  typename Selector< Edge >::value_type,
  std::ptrdiff_t,
  typename Selector< Edge >::reference,
  typename Selector< Edge >::pointer
  >
{
public:
  typedef Edge edge_type;
  typedef Selector< Edge > selector_type;
  typedef typename selector_type::value_type value_type;
  typedef typename selector_type::reference reference_type;
  typedef typename selector_type::pointer pointer_type;
  typedef typename selector_type::iterator underlying_iterator;
  typedef typename selector_type::ptr_type ptr_type;

  typedef PostOrder iterator;

private:
  ptr_type root_;
  bool at_top_;
  underlying_iterator pos_;
  std::deque< underlying_iterator > underlying_iterators_deque_;

private:
  PostOrder(ptr_type const& root, bool at_top);

public:
  static iterator begin(ptr_type const& root);
  static iterator end(ptr_type const& root);

public:
  ~PostOrder();

  reference_type operator *();
  pointer_type operator ->();

  iterator& operator ++();
  iterator operator ++(int);

private:
  void descend();

  template< typename E2 >
  friend bool operator ==(PostOrder< E2 > const& lhs, PostOrder< E2 > const& rhs);
};

template< typename Edge >
bool operator !=(PostOrder< Edge > const& lhs, PostOrder< Edge > const& rhs);

#include "iterator.hxx"

} // namespace iterator
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_ITERATOR_HPP_
