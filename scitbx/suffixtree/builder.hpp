#ifndef SUFFIXTREE_BUILDER_HPP_
#define SUFFIXTREE_BUILDER_HPP_

#include <scitbx/suffixtree/tree.hpp>

namespace scitbx
{

namespace suffixtree
{

namespace builder
{

template< typename Tree >
class Builder
{
public:
  typedef Tree tree_type;

  typedef typename tree_type::edge_type edge_type;
  typedef typename tree_type::edge_ptr_type edge_ptr_type;
  typedef typename tree_type::edge_weak_ptr_type edge_weak_ptr_type;

  typedef typename tree_type::word_ptr_type word_ptr_type;
  typedef typename tree_type::construction_ptr_type construction_ptr_type;

  typedef typename tree_type::word_type word_type;
  typedef typename tree_type::glyph_type glyph_type;
  typedef typename tree_type::index_type index_type;
  typedef typename tree_type::length_type length_type;
  typedef typename tree_type::word_length_ptr_type word_length_ptr_type;
  typedef typename tree_type::suffix_label_type suffix_label_type;

public:
  Builder();

protected:
  ~Builder();

protected:
  static edge_ptr_type const& tree_root(tree_type const& tree);
  static word_ptr_type const& tree_word_ptr(tree_type const& tree);
  static construction_ptr_type const& tree_construction_ptr(tree_type const& tree);
};

template< typename EdgePtr >
struct SuffixLinkerState
{
  typedef EdgePtr edge_ptr_type;
  virtual ~SuffixLinkerState();

  virtual bool process_existing(edge_ptr_type& stored, edge_ptr_type const& next) const= 0;
  virtual bool process_new(edge_ptr_type& stored, edge_ptr_type const& next) const = 0;
};

template< typename EdgePtr >
struct SuffixLinkerEmpty : public SuffixLinkerState< EdgePtr >
{
  typedef EdgePtr edge_ptr_type;

  virtual bool process_existing(edge_ptr_type& stored, edge_ptr_type const& next) const;
  virtual bool process_new(edge_ptr_type& stored, edge_ptr_type const& next) const;
};

template< typename EdgePtr >
struct SuffixLinkerPrimed : public SuffixLinkerState< EdgePtr >
{
  typedef EdgePtr edge_ptr_type;

  virtual bool process_existing(edge_ptr_type& stored, edge_ptr_type const& next) const;
  virtual bool process_new(edge_ptr_type& stored, edge_ptr_type const& next) const;
};

template< typename EdgePtr >
class SuffixLinker
{
public:
  typedef EdgePtr edge_ptr_type;
  typedef SuffixLinkerState< EdgePtr > state_type;
  typedef SuffixLinkerEmpty< EdgePtr > empty_state_type;
  typedef SuffixLinkerPrimed< EdgePtr > primed_state_type;

private:
  state_type const* state_;
  edge_ptr_type storage_;

  static empty_state_type const empty_state;
  static primed_state_type const primed_state;

public:
  SuffixLinker();
  ~SuffixLinker();

  void process_existing(const edge_ptr_type& next);
  void process_new(const edge_ptr_type& next);
};

template< typename Tree >
class Ukkonen : public Builder< Tree >
{
public:
  typedef Tree tree_type;
  typedef Builder< Tree > builder_type;

  typedef typename builder_type::edge_type edge_type;
  typedef typename builder_type::edge_ptr_type edge_ptr_type;
  typedef typename builder_type::edge_weak_ptr_type edge_weak_ptr_type;

  typedef typename builder_type::word_ptr_type word_ptr_type;
  typedef typename builder_type::construction_ptr_type construction_ptr_type;

  typedef typename builder_type::word_type word_type;
  typedef typename builder_type::glyph_type glyph_type;
  typedef typename builder_type::index_type index_type;
  typedef typename builder_type::length_type length_type;
  typedef typename builder_type::word_length_ptr_type word_length_ptr_type;
  typedef typename builder_type::suffix_label_type suffix_label_type;

  typedef SuffixLinker< edge_ptr_type > suffix_linker_type;
  typedef Cursor< edge_type, word_type > cursor_type;

private:
  edge_ptr_type tree_root_;
  word_ptr_type tree_word_ptr_;
  construction_ptr_type tree_construction_ptr_;

  cursor_type position_;
  suffix_label_type phase_;
  suffix_label_type extension_;
  bool is_attached_;

public:
  Ukkonen(tree_type const& tree);
  ~Ukkonen();

  bool is_attached() const;
  bool is_valid() const;

  void push_back(const glyph_type& glyph);
  void detach();
};

#include "builder.hxx"

} // namespace builder
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_BUILDER_HPP_
