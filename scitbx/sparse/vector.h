#ifndef SCITBX_SPARSE_VECTOR_H
#define SCITBX_SPARSE_VECTOR_H

#include <memory>
#include <algorithm>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <scitbx/error.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace sparse {


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
template<typename T>
class vector
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

  typedef af::shared<element> container_type;

  container_type elements;
  index_type size_;
  bool sorted;

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
    return std::binary_search(elements.begin(), elements.end(), element(i));
  }

  /// Number of non-zero elements
  index_type non_zeroes() const {
    compact();
    return elements.size();
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

  /// Fill the given dense vector.
  void fill_dense_vector(af::shared<T>& w) const {
    SCITBX_ASSERT(w.size() == size())
                 ( w.size() )( size() );
    for(const_iterator q =  begin(); q != end(); q++) {
      w[q.index()] = *q;
    }
  }

  /// The dense vector corresponding to this
  af::shared<T> as_dense_vector() const {
    af::shared<T> result(size(), 0.);
    fill_dense_vector(result);
    return result;
  }

  /// A copy of this, copying the elements.
  /** Since the copy constructor makes a shallow copy, thanks to
  the shallow copy semantic of af::shared, this is very necessary a
  member function. */
  vector deep_copy() const {
    vector result(size());
    result.elements = elements.deep_copy();
    result.sorted = sorted;
    return result;
  }

  /// Fill the given dense vector with a permutation of this
  template<class PermutationType>
  void fill_dense_vector_with_permutation(af::shared<T>& w,
                                          PermutationType const& perm) const
  {
    SCITBX_ASSERT(w.size() == perm.size() && perm.size() == size())
                 ( w.size() )( perm.size() )( size() );
    for(const_iterator q =  begin(); q != end(); q++) {
      w[ perm[q.index()] ] = *q;
    }
  }

  /// Whether this has been compacted
  bool is_compact() const { return sorted; }

  /// Perform summation and removal of duplicate indices, and sort indices.
  /** The record which was input last is kept in case of duplicate assignment.
  Return this object, for convenient chaining of operations. */
  vector const& compact() const {
    if (!sorted) const_cast<vector *>(this)->do_compact();
    return *this;
  }

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
};

}}

#endif
