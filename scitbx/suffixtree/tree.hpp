#ifndef SUFFIXTREE_TREE_HPP_
#define SUFFIXTREE_TREE_HPP_

#include <scitbx/suffixtree/edge.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>


#include <utility>
#include <iostream>

namespace scitbx
{

namespace suffixtree
{

// Forward declaration
namespace builder
{

template< typename Tree >
class Builder;

} // namespace builder

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
class Tree
{
public:
  typedef Word word_type;
  typedef boost::shared_ptr< word_type > word_ptr_type;
  typedef boost::shared_ptr< bool > construction_ptr_type;
  typedef typename word_type::glyph_type glyph_type;
  typedef typename word_type::index_type index_type;
  typedef typename word_type::length_type length_type;
  typedef typename word_type::const_length_ptr_type word_length_ptr_type;

  typedef SuffixLabel suffix_label_type;

  typedef edge::Edge<
    glyph_type,
    index_type,
    word_length_ptr_type,
    suffix_label_type,
    NodeAdaptor
    >
    edge_type;

  typedef typename edge_type::ptr_type edge_ptr_type;
  typedef typename edge_type::const_ptr_type const_edge_ptr_type;
  typedef typename edge_type::weak_ptr_type edge_weak_ptr_type;
  typedef typename edge_type::const_weak_ptr_type const_edge_weak_ptr_type;

private:
  edge_ptr_type root_;
  word_ptr_type word_ptr_;
  construction_ptr_type construction_ptr_;

public:
  Tree();
  ~Tree();

  const_edge_ptr_type root() const;
  word_type const& word() const;
  bool in_construction() const;

  friend class builder::Builder< Tree >;
};

template< typename Tree >
class Linkage
{
public:
  typedef Tree tree_type;
  typedef typename tree_type::edge_type edge_type;
  typedef typename tree_type::edge_ptr_type edge_ptr_type;
  typedef typename tree_type::glyph_type glyph_type;

public:
  static edge_ptr_type get_parent(edge_ptr_type const& edge_ptr);
  static edge_ptr_type get_suffix(edge_ptr_type const& edge_ptr);
  static edge_ptr_type get_child_with_label(
    edge_ptr_type const& edge_ptr,
    glyph_type const& glyph
    );
};

template< typename Tree >
class Movement
{
public:
  typedef Tree tree_type;

  typedef typename tree_type::edge_type edge_type;
  typedef typename tree_type::edge_ptr_type edge_ptr_type;
  typedef typename tree_type::word_type word_type;
  typedef typename tree_type::word_ptr_type word_ptr_type;
  typedef typename tree_type::glyph_type glyph_type;
  typedef typename tree_type::index_type index_type;
  typedef typename word_type::length_type length_type;
  typedef typename word_type::const_iterator word_iterator;
  typedef typename edge_type::iterator edge_iterator;

  typedef std::pair< edge_ptr_type, index_type > cursor_type;

  typedef Linkage< tree_type > linkage;

public:
  static edge_ptr_type const& get_edge_ptr(cursor_type const& cursor);
  static edge_ptr_type& get_edge_ptr(cursor_type& cursor);
  static index_type const& get_index(cursor_type const& cursor);
  static index_type& get_index(cursor_type& cursor);

  static bool is_at_edge_bottom(cursor_type const& cursor);

  static void set_to_edge_top(cursor_type& cursor, edge_ptr_type const& edge_ptr);
  static void forth(cursor_type& cursor);

  static cursor_type get_path_jump_destination(
    edge_ptr_type edge_ptr,
    word_iterator begin,
    word_iterator end
    );
  static cursor_type get_suffix_position(cursor_type const& cursor, word_type const& word);
};

#include "tree.hxx"

} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_TREE_HPP_
