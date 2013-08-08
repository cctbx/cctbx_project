#ifndef SUFFIXTREE_EDGE_HPP_
#define SUFFIXTREE_EDGE_HPP_

#include <scitbx/suffixtree/exception.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <utility>

namespace scitbx
{

namespace suffixtree
{

namespace edge
{

template< typename Key, typename Value >
struct ToConstPair
{
  typedef std::pair< const Key, boost::shared_ptr< Value > > const argument_type;
  typedef std::pair< const Key, boost::shared_ptr< Value const > > const result_type;

  result_type operator ()(argument_type const& arg) const;
};

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
  typedef boost::shared_ptr< edge_type const > const_ptr_type;

  typedef boost::weak_ptr< edge_type > weak_ptr_type;
  typedef boost::weak_ptr< edge_type const > const_weak_ptr_type;

  typedef Glyph glyph_type;
  typedef Index index_type;
  typedef WordLength word_length_type;
  typedef SuffixLabel suffix_label_type;
  typedef typename NodeAdapter< glyph_type, ptr_type >::type node_type;

  typedef ToConstPair< glyph_type, edge_type > converter_type;
  typedef typename node_type::iterator iterator;
  typedef typename node_type::value_type value_type;
  typedef boost::transform_iterator<
    converter_type,
    typename node_type::const_iterator
    >
    const_iterator;

public:
  Edge();
  virtual ~Edge();

  // Concrete functions
  bool empty() const;

  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  iterator find(const glyph_type& key);
  const_iterator find(const glyph_type& key) const;

  std::pair< iterator, bool > insert(value_type const& val);

  // Virtual functions
  virtual index_type const& start() const = 0;
  virtual index_type& start() = 0;

  virtual index_type stop() const = 0;

  virtual suffix_label_type const& label() const = 0;

  virtual const_weak_ptr_type parent() const = 0;
  virtual weak_ptr_type& parent() = 0;

  virtual const_weak_ptr_type suffix() const = 0;
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
  typedef typename edge_type::const_ptr_type const_ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;
  typedef typename edge_type::const_weak_ptr_type const_weak_ptr_type;

private:
  node_type node_;

private:
  static index_type const shared_start;

public:
  Root();
  virtual ~Root();

  virtual index_type const& start() const;
  virtual index_type& start();

  virtual index_type stop() const;

  virtual suffix_label_type const& label() const;

  virtual const_weak_ptr_type parent() const;
  virtual weak_ptr_type& parent();

  virtual const_weak_ptr_type suffix() const;
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
  typedef typename edge_type::const_ptr_type const_ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;
  typedef typename edge_type::const_weak_ptr_type const_weak_ptr_type;

private:
  index_type start_;
  index_type stop_;
  node_type node_;
  weak_ptr_type parent_ptr_;
  weak_ptr_type suffix_ptr_;

public:
  Branch(index_type const& start, index_type const& end);
  virtual ~Branch();

  virtual index_type const& start() const;
  virtual index_type& start();

  virtual index_type stop() const;

  virtual suffix_label_type const& label() const;

  virtual const_weak_ptr_type parent() const;
  virtual weak_ptr_type& parent();

  virtual const_weak_ptr_type suffix() const;
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
  typedef typename edge_type::const_ptr_type const_ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;
  typedef typename edge_type::const_weak_ptr_type const_weak_ptr_type;

private:
  index_type start_;
  word_length_type word_length_;
  suffix_label_type suffix_label_;
  weak_ptr_type parent_ptr_;

  static node_type const shared_node;

public:
  Leaf(
    index_type const& start,
    word_length_type const& word_length,
    suffix_label_type const& suffix_label
    );
  virtual ~Leaf();

  virtual index_type const& start() const;
  virtual index_type& start();

  virtual index_type stop() const;

  virtual suffix_label_type const& label() const;

  virtual const_weak_ptr_type parent() const;
  virtual weak_ptr_type& parent();

  virtual const_weak_ptr_type suffix() const;
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
