#ifndef SUFFIXTREE_WORD_HPP_
#define SUFFIXTREE_WORD_HPP_

#include <boost/range/iterator_range.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <iterator>
#include <vector>

namespace scitbx
{

namespace suffixtree
{

namespace word
{

template< typename Glyph >
class Single
{
public:
  typedef Glyph glyph_type;
  typedef std::vector< glyph_type > data_type;
  typedef typename data_type::const_iterator const_iterator;
  typedef typename data_type::size_type index_type;
  typedef typename data_type::size_type length_type;
  typedef length_type const const_length_type;
  typedef boost::shared_ptr< length_type > length_ptr_type;
  typedef boost::shared_ptr< const_length_type > const_length_ptr_type;
  typedef boost::iterator_range< const_iterator > substring_type;
  typedef Single type;

private:
  data_type data_;
  length_ptr_type length_ptr_;

public:
  Single();
  ~Single();

  inline void push_back(glyph_type const& ch);

  const_length_ptr_type length_ptr() const;
  length_type size() const;
  substring_type substring(index_type const& begin, index_type const& end) const;

  const_iterator get_iterator_to(index_type const& index) const;

  glyph_type const& operator [](index_type const& index) const;
  glyph_type& operator [](index_type const& index);
};

#include "word.hxx"

} // namespace word
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_WORD_HPP_
