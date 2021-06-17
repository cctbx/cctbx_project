#ifndef SUFFIXTREE_TREE_HPP_
#define SUFFIXTREE_TREE_HPP_

#include <scitbx/suffixtree/edge.hpp>
#include <scitbx/suffixtree/exception.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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

template< typename Edge, typename Word >
class Cursor
{
public:
  typedef Edge edge_type;
  typedef typename edge_type::ptr_type edge_ptr_type;
  typedef typename edge_type::weak_ptr_type edge_weak_ptr_type;
  typedef typename edge_type::iterator edge_iterator;
  typedef typename edge_type::glyph_type glyph_type;
  typedef typename edge_type::index_type index_type;

  typedef Word word_type;
  typedef boost::shared_ptr< word_type const > word_ptr_type;
  typedef typename word_type::length_type length_type;
  typedef typename word_type::const_iterator word_iterator;

private:
  word_ptr_type word_ptr_;
  edge_ptr_type edge_ptr_;
  index_type index_;

public:
  Cursor(edge_ptr_type const& edge_ptr, word_ptr_type const& word_ptr);
  ~Cursor();

  edge_ptr_type const& get_edge_ptr() const;
  index_type const& get_index() const;

  bool is_at_edge_top() const;
  bool is_at_edge_bottom() const;
  bool is_at_leaf_bottom() const;

  glyph_type const& get_current_character() const;
  word_type const& get_word() const;
  //edge_ptr_type get_child_with_label(glyph_type const& glyph) const;

  void forth_with(glyph_type const& glyph);
  void break_edge_here();
  void to_suffix_position();

private:
  void forth_to_child(glyph_type const& label);
  void forth_on_edge();
  void path_jump_from_top_of(
    edge_ptr_type edge_ptr,
    word_iterator begin,
    word_iterator end
    );
};

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
  typedef typename edge_type::weak_ptr_type edge_weak_ptr_type;

  typedef Cursor< edge_type, word_type > cursor_type;

private:
  edge_ptr_type root_;
  word_ptr_type word_ptr_;
  construction_ptr_type construction_ptr_;

public:
  Tree();
  ~Tree();

  edge_ptr_type root() const;
  word_type const& word() const;
  bool in_construction() const;
  cursor_type cursor() const;

  friend class builder::Builder< Tree >;

  template< typename W, typename SL, template< typename, typename > class NA >
  friend bool operator ==(Tree< W, SL, NA > const& lhs, Tree< W, SL, NA > const& rhs);
};

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
bool operator ==(
  Tree< Word, SuffixLabel, NodeAdaptor > const& lhs,
  Tree< Word, SuffixLabel, NodeAdaptor > const& rhs
  );

#include "tree.hxx"

} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_TREE_HPP_
