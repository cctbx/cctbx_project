#ifndef SUFFIXTREE_EDGE_HPP_
#define SUFFIXTREE_EDGE_HPP_

#include <scitbx/suffixtree/exception.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include <utility>
#include <iostream>

namespace scitbx
{

namespace suffixtree
{

namespace edge
{

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
class Edge
{
public:
  typedef Edge edge_type;
  typedef boost::shared_ptr< edge_type > ptr_type;
  typedef boost::weak_ptr< edge_type > weak_ptr_type;

  typedef Glyph glyph_type;
  typedef Index index_type;
  typedef WordLength word_length_type;
  typedef SuffixLabel suffix_label_type;
  typedef typename NodeAdapter< glyph_type, ptr_type >::type node_type;

  typedef typename node_type::iterator iterator;
  typedef typename node_type::value_type value_type;
  typedef typename node_type::const_iterator const_iterator;

public:
  Edge();
  virtual ~Edge();

  // Concrete functions
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  iterator find(const glyph_type& key);
  const_iterator find(const glyph_type& key) const;

  std::pair< iterator, bool > insert(value_type const& val);

  // Convenience functions
  ptr_type get_parent() const;
  ptr_type get_suffix() const;

  void set_parent(ptr_type const& parent);

private:
  ptr_type promote_link(weak_ptr_type const& ptr) const;

public:
  ptr_type get_child_with_label(glyph_type const& label) const;

  void attach_child(ptr_type const& child, glyph_type const& label);
  bool attach_child_if_not_present(ptr_type const& child, glyph_type const& label);

  // Virtual functions
  virtual bool empty() const = 0;

  virtual index_type get_start() const = 0;
  virtual void set_start(index_type const& value) = 0;

  virtual index_type get_stop() const = 0;

  virtual suffix_label_type const& label() const = 0;

  virtual weak_ptr_type const& parent() const = 0;
  virtual weak_ptr_type& parent() = 0;

  virtual weak_ptr_type const& suffix() const = 0;
  virtual weak_ptr_type& suffix() = 0;

  virtual bool is_root() const = 0;
  virtual bool is_leaf() const = 0;

private:
  virtual node_type const& node() const = 0;
  virtual node_type& node() = 0;

public:
  static ptr_type root();
  static ptr_type branch(index_type const& start, index_type const& end);
  static ptr_type leaf(
    index_type const& start,
    word_length_type const& word_length,
    suffix_label_type const& suffix_label
    );
};

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
class Root : public Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >
{
public:
  typedef Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > edge_type;

  typedef typename edge_type::glyph_type glyph_type;
  typedef typename edge_type::index_type index_type;
  typedef typename edge_type::word_length_type word_length_type;
  typedef typename edge_type::suffix_label_type suffix_label_type;
  typedef typename edge_type::node_type node_type;
  typedef typename edge_type::ptr_type ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;

private:
  node_type node_;

public:
  Root();
  virtual ~Root();

  virtual bool empty() const;

  virtual index_type get_start() const;
  virtual void set_start(index_type const& value);

  virtual index_type get_stop() const;

  virtual suffix_label_type const& label() const;

  virtual weak_ptr_type const& parent() const;
  virtual weak_ptr_type& parent();

  virtual weak_ptr_type const& suffix() const;
  virtual weak_ptr_type& suffix();

  virtual bool is_root() const;
  virtual bool is_leaf() const;

private:
  virtual node_type const& node() const;
  virtual node_type& node();
};

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
class Branch : public Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >
{
public:
  typedef Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > edge_type;

  typedef typename edge_type::glyph_type glyph_type;
  typedef typename edge_type::index_type index_type;
  typedef typename edge_type::word_length_type word_length_type;
  typedef typename edge_type::suffix_label_type suffix_label_type;
  typedef typename edge_type::node_type node_type;
  typedef typename edge_type::ptr_type ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;

private:
  index_type start_;
  index_type stop_;
  node_type node_;
  weak_ptr_type parent_ptr_;
  weak_ptr_type suffix_ptr_;

public:
  Branch(index_type const& start, index_type const& end);
  virtual ~Branch();

  virtual bool empty() const;

  virtual index_type get_start() const;
  virtual void set_start(index_type const& value);

  virtual index_type get_stop() const;

  virtual suffix_label_type const& label() const;

  virtual weak_ptr_type const& parent() const;
  virtual weak_ptr_type& parent();

  virtual weak_ptr_type const& suffix() const;
  virtual weak_ptr_type& suffix();

  virtual bool is_root() const;
  virtual bool is_leaf() const;

private:
  virtual node_type const& node() const;
  virtual node_type& node();
};

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
class Leaf : public Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter >
{
public:
  typedef Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > edge_type;

  typedef typename edge_type::glyph_type glyph_type;
  typedef typename edge_type::index_type index_type;
  typedef typename edge_type::word_length_type word_length_type;
  typedef typename edge_type::suffix_label_type suffix_label_type;
  typedef typename edge_type::node_type node_type;
  typedef typename edge_type::ptr_type ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;

private:
  index_type start_;
  word_length_type word_length_;
  suffix_label_type suffix_label_;
  weak_ptr_type parent_ptr_;

public:
  Leaf(
    index_type const& start,
    word_length_type const& word_length,
    suffix_label_type const& suffix_label
    );
  virtual ~Leaf();

  virtual bool empty() const;

  virtual index_type get_start() const;
  virtual void set_start(index_type const& value);

  virtual index_type get_stop() const;

  virtual suffix_label_type const& label() const;

  virtual weak_ptr_type const& parent() const;
  virtual weak_ptr_type& parent();

  virtual weak_ptr_type const& suffix() const;
  virtual weak_ptr_type& suffix();

  virtual bool is_root() const;
  virtual bool is_leaf() const;

private:
  virtual node_type const& node() const;
  virtual node_type& node();
};

#include "edge.hxx"

} // namespace edge
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_EDGE_HPP_
