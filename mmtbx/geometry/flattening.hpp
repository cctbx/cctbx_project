// Based on http://stackoverflow.com/questions/3623082/flattening-iterator
#ifndef MMTBX_GEOMETRY_FLATTEN_H
#define MMTBX_GEOMETRY_FLATTEN_H

#include <iterator>
#include <vector>

// A forward iterator that "flattens" a container of containers.  For example,
// a vector<vector<int>> containing { { 1, 2, 3 }, { 4, 5, 6 } } is iterated as
// a single range, { 1, 2, 3, 4, 5, 6 }.
namespace mmtbx
{

namespace geometry
{

namespace utility
{

template< typename OuterIterator, typename InnerIterator >
class flattening_iterator
{
public:

    typedef OuterIterator  outer_iterator;
    typedef InnerIterator  inner_iterator;

    typedef std::forward_iterator_tag                iterator_category;
    typedef typename inner_iterator::value_type      value_type;
    typedef typename inner_iterator::difference_type difference_type;
    typedef typename inner_iterator::pointer         pointer;
    typedef typename inner_iterator::reference       reference;

    flattening_iterator() { }
    flattening_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
    flattening_iterator(outer_iterator it, outer_iterator end)
        : outer_it_(it),
          outer_end_(end)
    {
        if (outer_it_ == outer_end_) { return; }

        inner_it_ = outer_it_->begin();
        advance_past_empty_inner_containers();
    }

    reference operator*()  const { return *inner_it_;  }
    pointer   operator->() const { return &*inner_it_; }

    flattening_iterator& operator++()
    {
        ++inner_it_;
        if (inner_it_ == outer_it_->end())
            advance_past_empty_inner_containers();
        return *this;
    }

    flattening_iterator operator++(int)
    {
        flattening_iterator it(*this);
        ++*this;
        return it;
    }

    friend bool operator==(const flattening_iterator& a,
                           const flattening_iterator& b)
    {
        if (a.outer_it_ != b.outer_it_)
            return false;

        if (a.outer_it_ != a.outer_end_ &&
            b.outer_it_ != b.outer_end_ &&
            a.inner_it_ != b.inner_it_)
            return false;

        return true;
    }

    friend bool operator!=(const flattening_iterator& a,
                           const flattening_iterator& b)
    {
        return !(a == b);
    }

private:

    void advance_past_empty_inner_containers()
    {
        while (outer_it_ != outer_end_ && inner_it_ == outer_it_->end())
        {
            ++outer_it_;
            if (outer_it_ != outer_end_)
                inner_it_ = outer_it_->begin();
        }
    }

    outer_iterator outer_it_;
    outer_iterator outer_end_;
    inner_iterator inner_it_;
};

// Flattened range
template< typename Range >
struct flattening_range
{
public:
  typedef Range range_type;
  typedef typename range_type::value_type value_type;
  typedef std::vector< range_type > storage_type;
  typedef flattening_iterator<
    typename storage_type::iterator,
    typename storage_type::iterator::value_type::iterator
    >
    iterator;
  typedef flattening_iterator<
    typename storage_type::const_iterator,
    typename storage_type::const_iterator::value_type::const_iterator
    >
    const_iterator;

  storage_type storage;

  iterator begin()
  {
    return iterator( storage.begin(), storage.end() );
  }

  iterator end()
  {
    return iterator( storage.end(), storage.end() );
  }

  const_iterator begin() const
  {
    return const_iterator( storage.begin(), storage.end() );
  }

  const_iterator end() const
  {
    return const_iterator( storage.end(), storage.end() );
  }
};

} // namespace utility
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_FLATTEN_H
