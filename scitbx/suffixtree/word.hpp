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
  typedef std::size_t index_type;
  typedef std::size_t length_type;
  typedef boost::shared_ptr< length_type > length_ptr_type;
  typedef boost::shared_ptr< const index_type > const_length_ptr_type;
  typedef boost::iterator_range< const_iterator > substring_type;
  typedef Single type;

private:
  data_type data_;
  length_ptr_type length_ptr_;

public:
  Single();
  ~Single();

  inline void push_back(const glyph_type& ch);

  const_length_ptr_type length_ptr() const;
  index_type size() const;
  substring_type substring(const index_type& begin, const index_type& end) const;

  const_iterator get_iterator_to(const index_type& index) const;

  const glyph_type& operator [](const index_type& index) const;
  glyph_type& operator [](const index_type& index);
};

template< typename Glyph, typename Traits >
class Multiple
{
public:
  typedef Glyph glyph_type;
  typedef Traits traits_type;
  typedef std::vector< glyph_type > data_type;
  typedef typename data_type::const_iterator const_iterator;
  typedef std::size_t index_type;
  typedef boost::shared_ptr< index_type > length_ptr_type;
  typedef boost::shared_ptr< const index_type > length_type;
  typedef boost::iterator_range< const_iterator > substring_type;

private:
  data_type data_;
  length_ptr_type length_ptr_;

public:
  Multiple();
  ~Multiple();

  inline void push_back(const glyph_type& glyph);

  length_type length() const;
  substring_type substring(const index_type& begin, const index_type& end) const;

  const glyph_type& operator [](const index_type& index) const;
  glyph_type& operator [](const index_type& index);

private:
  void set_word_boundary();
};

#include "word.hxx"

} // namespace word
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_WORD_HPP_
