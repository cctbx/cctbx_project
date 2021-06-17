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
struct IteratorTraits
{
  typedef Edge edge_type;
  typedef typename edge_type::const_iterator underlying_const_iterator;
  typedef typename std::forward_iterator_tag iterator_category;
  typedef typename edge_type::ptr_type ptr_type;
  typedef ptr_type const value_type;

  typedef std::ptrdiff_t difference_type;

  typedef value_type& reference;
  typedef value_type* pointer;
};

template< typename Edge >
class PreOrder
{
public:
  typedef Edge edge_type;
  typedef IteratorTraits< Edge > traits_type;
  typedef typename traits_type::iterator_category iterator_category;
  typedef typename traits_type::value_type value_type;
  typedef typename traits_type::difference_type difference_type;
  typedef typename traits_type::reference reference;
  typedef typename traits_type::pointer pointer;
  typedef typename traits_type::underlying_const_iterator
    underlying_const_iterator;
  typedef typename traits_type::ptr_type ptr_type;

  typedef PreOrder iterator;

private:
  ptr_type root_;
  bool at_top_;
  underlying_const_iterator pos_;
  std::deque< underlying_const_iterator > underlying_iterators_deque_;

private:
  PreOrder(ptr_type const& root, bool at_top);

public:
  static iterator begin(ptr_type const& root);
  static iterator end(ptr_type const& root);

public:
  ~PreOrder();

  reference operator *() const;
  pointer operator ->() const;

  iterator& operator ++();
  iterator operator ++(int);

  template< typename E2 >
  friend bool operator ==(PreOrder< E2 > const& lhs, PreOrder< E2 > const& rhs);
};

template< typename Edge >
bool operator !=(PreOrder< Edge > const& lhs, PreOrder< Edge > const& rhs);

template< typename Edge >
class PostOrder
{
public:
  typedef Edge edge_type;
  typedef IteratorTraits< Edge > traits_type;
  typedef typename traits_type::iterator_category iterator_category;
  typedef typename traits_type::value_type value_type;
  typedef typename traits_type::difference_type difference_type;
  typedef typename traits_type::reference reference;
  typedef typename traits_type::pointer pointer;
  typedef typename traits_type::underlying_const_iterator
    underlying_const_iterator;
  typedef typename traits_type::ptr_type ptr_type;

  typedef PostOrder iterator;

private:
  ptr_type root_;
  bool at_top_;
  underlying_const_iterator pos_;
  std::deque< underlying_const_iterator > underlying_iterators_deque_;

private:
  PostOrder(ptr_type const& root, bool at_top);

public:
  static iterator begin(ptr_type const& root);
  static iterator end(ptr_type const& root);

public:
  ~PostOrder();

  reference operator *() const;
  pointer operator ->() const;

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
