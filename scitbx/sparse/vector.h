#ifndef SCITBX_SPARSE_VECTOR_H
#define SCITBX_SPARSE_VECTOR_H

#include <memory>
#include <algorithm>
#include <functional>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <scitbx/error.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace sparse {


template <class E>
struct vector_expression
{
  E const &operator()() const {return static_cast<E const &>(*this); }
};


/** A sparse vector represented as a sequence of records,
 each containing the value and the index of a non-zero element.

 The semantic is as follow for a vector v.

 (1) If no value has been assigned to v[i], then v[i] == 0 and no data is stored
 corresponding to that index i.

 (2) After an assignment v[i] = x, even if x is zero, a pair (i,x) is stored and
 v[i] == x

 In sparse algorithm, v[i] == 0 in case (1) corresponds to structural zeroes:
 those elements are never touched by the algorithm; whereas v[i] == 0 in case (2)
 results from the assignement to v[i] of an expression which happens to be zero:
 that's a coincidential cancellation.

 Successive assignments and augmented assignments work as expected, i.e.

 @code
 sparse::vector v(3);
 v[i] += 1 // v[i] == 1
 ...
 v[i] += 2 // v[i] == 3
 ...
 v[i] = 4  // v[i] == 4
 ...
 v[i] -= 1 // v[i] == 3
 ...
 v[i] = 6  // v[i] == 6
 @endcode

 Such a sequence of assignments never fetches the value v[i]:
 v records the values to assign, add or substract as pairs (index, value)
 in the order they come. This efficiency is only available in C++:
 in Python, v[i] += ... will fetch the value v[i].

 (3) Many operations require that elements are sorted by increasing
 index without duplicate indices, a layout we will referred to as "compact" in
 the following. This layout can be achieved by calling compact() beforehand
 and it is automatically called by default when needed:
 one of the most important example being

 @code double x = v[i] @endcode

 (4) The pre-condition that v[i] = ... is only possible if i is less than size()
 is not enforced for efficiency reasons. Calling compact() will however prune
 those illegal elements.

 Implementation note:
 The C++ standard rules that the private types and members of a class are not
 accessible to its nested classes. A defect report (issue 45, [1]) has however
 overturned that decision. It does not have the status "TC1" and is therefore
 not part of the C++ standard yet. However, all compilers tested with the cctbx
 at the time of writing have implemented the recommendation of issue 45,
 i.e. that nested classes have access to all members of the outer class.
 It should be noted that cxx on Tru64 is true to the standard and does not
 implement it but this platform is not longer tested against.

[1] http://www.open-std.org/jtc1/sc22/wg21/docs/cwg_defects.html
*/
template<typename T, template<class> class ContainerType=af::shared>
class vector : public vector_expression< vector<T, ContainerType> >,
               public af::expression< vector<T, ContainerType> >,
               public boost::equality_comparable< vector<T, ContainerType> >
{
public:
  typedef T value_type;
  typedef std::size_t index_type;
  typedef std::ptrdiff_t index_difference_type;
  typedef af::const_ref<value_type> dense_vector_const_ref;

private:
  /// An element (index, value) of the sparse vector
  /** The highest bit of index is used to record whether the value is
   to be assigned or added.
   */
  class element : boost::totally_ordered<element>
  {
  private:
    index_type index_;
    value_type value_;

    static
    index_type sum_flag()
    {
      return index_type(1u) << (8*sizeof(index_type) - 1);
    }

  public:
    /// Construct an element to assign
    element(index_type i, value_type x=0)
        : index_(i & ~sum_flag()), value_(x)
    {}

    /// Construct an element to add to a sum
    element(index_type i, bool sum, value_type x=0)
        : index_(i | sum_flag()), value_(x)
    {}

    index_type index() const {
      return index_ & ~sum_flag();
    }

    value_type value() const { return value_; }

    value_type &value() { return value_; }

    bool summed() { return index_ & sum_flag(); }

    template <class PermutationType>
    void apply(PermutationType const &p) {
      index_ = p[index()] | (index_ & sum_flag());
    }

    bool
    operator==(element const &other) const { return index() == other.index(); }

    bool
    operator<(element const &other) const { return index() < other.index(); }
  };

  typedef ContainerType<element> container_type;

  container_type elements;
  bool sorted;
  index_type size_;

  value_type get(index_type i) const {
    compact();
    if (is_structurally_zero()) return 0;
    typename container_type::const_iterator
    p = std::lower_bound(elements.begin(), elements.end(), element(i));
    if (p != elements.end() && p->index() == i) return p->value();
    else return 0;
  }

  void set(index_type i, value_type x) {
    elements.push_back(element(i, x));
    sorted = false;
  }

  void add(index_type i, value_type x) {
    elements.push_back(element(i, true, x));
    sorted = false;
  }

  void do_compact() {
    if (elements.size()) {
      std::stable_sort(elements.begin(), elements.end());

      typedef typename container_type::iterator iter_t;
      iter_t q = elements.end()-1, overwrite = q;
      while (q >= elements.begin())
      {
        iter_t p;
        index_type current_index = q->index();
        if (current_index >= size_) {
          --q;
          continue;
        }
        for (p=q; p >= elements.begin() + 1; --p) {
          if (p[-1].index() != current_index) break;
          if (!p->summed()) break;
        }
        value_type x = p->value();
        for (iter_t r=p+1; r <= q; ++r) x += r->value();
        *overwrite-- = element(current_index, x);
        for (q=p-1; q >= elements.begin(); --q) {
          if (q->index() != current_index) break;
        }
      }
      elements.erase(elements.begin(), overwrite + 1);
    }
    sorted = true;
  }

public:
  /// Const iterator over the records
  class const_iterator
    : public boost::forward_iterator_helper<const_iterator, value_type const>
  {
    typename container_type::const_iterator p;

  public:
    const_iterator()
    {}

    const_iterator(typename container_type::const_iterator q)
      : p(q)
    {}

    value_type operator*() const {
      return p->value();
    }

    bool operator==(const_iterator const& i) const {
      return p == i.p;
    }

    const_iterator& operator++() {
      p++;
      return *this;
    }

    index_type index() const {
      return p->index();
    }
  };

  /// Iterator over the records
  class iterator
    : public boost::forward_iterator_helper<iterator, value_type>
  {
    typename container_type::iterator p;

  public:
    iterator()
    {}

    iterator(typename container_type::iterator q)
      : p(q)
    {}

    value_type& operator*() const {
      return p->value();
    }

    bool operator==(iterator const& i) const {
      return p == i.p;
    }

    iterator& operator++() {
      p++;
      return *this;
    }

    index_type index() const {
      return p->index();
    }
  };

  /// A const reference to an element of given index
  /** This the type of object returned by v[i] for a const vector v
      and an index i
  */
  class element_const_reference
  {
  private:
    vector const &v;
    index_type i;

  public:
    /// Construct the reference to the element of index j in u
    element_const_reference(vector const &u, index_type j) : v(u), i(j)
    {}

    /* Without the destructor, the constructor would not be inlined
    by g++ 4.0 */
    ~element_const_reference()
    {}

    /// Triggered by using v[i] in an expression
    operator value_type() const {
      return v.get(i);
    }
  };


  /// A reference to an element of given index
  /** This the type of object returned by v[i] for a vector v and an index i
  */
  class element_reference
  {
  private:
    vector &v;
    index_type i;

  public:
    /// Construct the reference to the element of index j in u
    element_reference(vector &u, index_type j) : v(u), i(j)
    {}

    /* Without the destructor, the constructor would not be inlined
    by g++ 4.0 */
    ~element_reference()
    {}

    /// Triggered by using v[i] in an expression
    /** Runtime scales as O(Log n) where n is the number of non-zero elements,
        if the vector is compacted. Otherwise, O(Log n) plus the cost of
        compacting it.
     */
    operator value_type() const {
      return v.get(i);
    }

    /// Triggered by an assignment v[i] = ...
    /** Runtime scales as O(1)
     */
    element_reference &operator=(value_type x) {
      v.set(i, x);
      return *this;
    }

    /// Triggered by an assignment v[i] += ...
    /** Runtime scales as O(1)
     */
    element_reference &operator+=(value_type x) {
      v.add(i, x);
      return *this;
    }

    /// Triggered by an assignment v[i] -= ...
    /** Runtime scales as O(1)
     */
    element_reference &operator-=(value_type x) {
      v.add(i, -x);
      return *this;
    }
  };

  /// Construct a zero vector of size n
  vector(index_type n) : sorted(false), size_(n) {}

  /// Construct from a generic expression
  template <class E>
  vector(vector_expression<E> const &e)
    : size_(e().size())
  {
    e().assign_to(*this);
  }

  /// Assignment of a generic expression
  template <class E>
  vector &operator=(vector_expression<E> const &e) {
    size_ = e().size();
    elements.clear();
    e().assign_to(*this);
    return *this;
  }

  /// Assign to the given reference to a dense vector.
  /** Along with member function size() this allows
   \begincode
   af::shared<T> u = v;
   af::ref<T> u = v;
   \endcode
   for any sparse::vector<T> v
   */
  void assign_to(af::ref<T> const &w) const {
    SCITBX_ASSERT(w.size() == size())
    ( w.size() )( size() );
    for(const_iterator q =  begin(); q != end(); q++) {
      w[q.index()] = *q;
    }
  }

  /// Specify whether the vector shall be considered compacted or not
  /** set_compact(true) is very dangerous: it shall only be used
      along with algorithms provingly building vectors with increasing
      unique indices.
   */
  void set_compact(bool f) const { const_cast<vector &>(*this).sorted = true; }

  /// Whether this has been compacted
  bool is_compact() const { return sorted; }

  /// Perform summation and removal of duplicate indices, and sort indices.
  /** The record which was input last is kept in case of duplicate assignment.
  Return this object, for convenient chaining of operations. */
  vector const& compact() const {
    if (!sorted) const_cast<vector *>(this)->do_compact();
    return *this;
  }

  /// An iterator pointing to the first record
  const_iterator begin() const {
    return const_iterator(elements.begin());
  }

  /// An iterator pointing past the last record
  const_iterator end() const {
    return const_iterator(elements.end());
  }

  /// An iterator pointing to the first record
  iterator begin() {
    return iterator(elements.begin());
  }

  /// An iterator pointing past the last record
  iterator end() {
    return iterator(elements.end());
  }

  /// Dimension of the vector, i.e. number of zero or non-zero elements
  index_type size() const {
    return size_;
  }

  /// Whether there is no potential non-zero elements
  bool is_structurally_zero() const {
    return elements.size() == 0;
  }

  /// Whether the element of index i is a structural zero
  bool is_structural_zero(index_type i) const {
    compact();
    return !std::binary_search(elements.begin(), elements.end(), element(i));
  }

  /// Equality
  bool operator==(vector const &other) const {
    compact();
    other.compact();
    return elements.size() == other.elements.size()
        && std::equal(elements.begin(), elements.end(), other.elements.begin());
  }

  /// Number of non-zero elements
  index_type non_zeroes() const {
    compact();
    return elements.size();
  }

  /// Zero vector
  void zero() {
    elements.clear();
    sorted = true;
  }

  /// Subscripting
  element_const_reference operator[](index_type i) const {
    return element_const_reference(*this, i);
  }

  /// Subscripting
  /** Assignment v[i] = ... may introduce a duplicate index, a problem
  solved for the getter v[i] by calling compact() before searching
  that index.
  */
  element_reference operator[](index_type i) {
    return element_reference(*this, i);
  }

  /// Dense row vector times sparse vector
  friend value_type operator*(dense_vector_const_ref const &u,
                              vector const &v)
  {
    v.compact();
    value_type result = 0;
    for (const_iterator pv=v.begin(); pv != v.end(); ++pv) {
      index_type i = pv.index();
      value_type u_i = u[i];
      value_type v_i = *pv;
      result += u_i * v_i;
    }
    return result;
  }

  /// sparse vector times dense column vector
  friend value_type operator*(vector const &u,
                              dense_vector_const_ref const &v)
  {
    return v*u;
  }

  /// The dense vector corresponding to this
  af::shared<T> as_dense_vector() const {
    af::shared<T> result(size(), 0.);
    result.ref() = *this;
    return result;
  }

  /// Clone of this, copying the element
  /** Even if the ContainerType used to store elements has a sharing semantic.
   */
  vector clone() const {
    vector v(size());
    BOOST_FOREACH(element const &e, elements) v.elements.push_back(e);
    v.sorted = sorted;
    return v;
  }

  /// Linear algebra
  //@{
  friend value_type operator*(vector const &u, vector const &v) {
    return dot_product(u, v).result;
  }

  friend
  vector operator+(vector const &u, vector const &v) {
    return vector_op_vector< std::plus<T> >(u, v).result;
  }

  friend
  vector operator-(vector const &u, vector const &v) {
    return vector_op_vector< std::minus<T> >(u, v).result;
  }

  vector operator-() const {
    vector v(size());
    for (const_iterator p=begin(); p != end(); ++p) {
      v[p.index()] = -(*p);
    }
    return v;
  }

  vector &operator*=(T a) {
    BOOST_FOREACH(element &e, elements) e.value() *= a;
    return *this;
  }

  vector &operator/=(T a) {
    return (*this) *= 1./a;
  }

  friend
  vector operator*(vector const &u, T a) {
    vector v = u.clone();
    v *= a;
    return v;
  }

  friend
  vector operator*(T a, vector const &u) { return u*a; }

  friend
  vector operator/(vector const &u, T a) {
    vector v = u.clone();
    v /= a;
    return v;
  }
  //@}

  /// Permute the elements of this, in place
  /** Return this object, for convenient chaining of operations */
  template<class PermutationType>
  vector& permute(PermutationType const& permutation) {
    SCITBX_ASSERT(size() == permutation.size())
                 ( size() )( permutation.size() );
    BOOST_FOREACH(element &e, elements) {
      e.apply(permutation);
    }
    return *this;
  }

private:
  template <class Heir>
  struct vector_op_vector_core
  {
    void loop(vector const &u, vector const &v)
    {
      SCITBX_ASSERT(u.size() == v.size())( u.size() )( v.size() );
      u.compact();
      v.compact();
      Heir &self = static_cast<Heir &>(*this);
      for(const_iterator p=u.begin(), q=v.begin();;) {
        if (p == u.end()) {
          for(; q != v.end(); ++q) self(q.index(), 0, *q);
          break;
        }
        else if (q == v.end()) {
          for(; p != u.end(); ++p) self(p.index(), *p, 0);
          break;
        }
        std::size_t i=p.index(), j=q.index();
        if      (i < j) self(i, *p++, 0);
        else if (i > j) self(j, 0, *q++);
        else            self(i, *p++, *q++);
      }
    }
  };

  template <class OperatorType>
  struct vector_op_vector
  : vector_op_vector_core<vector_op_vector< OperatorType> >
  {
    vector result;
    OperatorType op;

    vector_op_vector(vector const &u, vector const &v)
    : result(u.size())
    {
      this->loop(u, v);
      result.set_compact(true); // by construction
    }

    void operator()(index_type i, value_type x, value_type y) {
      result[i] = op(x, y);
    }
  };

  struct dot_product : vector_op_vector_core<dot_product>
  {
    value_type result;

    dot_product(vector const &u, vector const &v)
    : result(0)
    {
      this->loop(u,v);
    }

    void operator()(index_type i, value_type x, value_type y) {
      result += x*y;
    }
  };

  template <typename VectorType, class PermutationType>
  friend struct permuted;
};


/// Permuted sparse vector expression
template <typename VectorType, class PermutationType>
struct permuted :
  af::expression< permuted<VectorType, PermutationType> >,
  vector_expression< permuted<VectorType, PermutationType> >
{
  VectorType const &v;
  PermutationType const &permutation;

  permuted(VectorType const &v, PermutationType const &p)
  : v(v), permutation(p)
  {
    SCITBX_ASSERT(v.size() == p.size())( v.size() )( p.size() );
  }

  std::size_t size() const { return v.size(); }

  /// Assign to a dense vector
  void assign_to(af::ref<typename VectorType::value_type> const &w) const {
    SCITBX_ASSERT(w.size() == v.size())( w.size() )( v.size() );
    for(typename VectorType::const_iterator q = v.begin(); q != v.end(); q++) {
      w[ permutation[q.index()] ] = *q;
    }
  }

  /// Assign to a sparse vector
  void assign_to(VectorType &w) const {
    w = v;
    w.permute(permutation);
  }
};

/// Permutation of this as an expression
template <typename VectorType, class PermutationType>
permuted<VectorType, PermutationType>
permute(VectorType const &v, PermutationType const &p) {
  return permuted<VectorType, PermutationType>(v, p);
}

}}

#endif
