#ifndef SCITBX_SPARSE_VECTOR_H
#define SCITBX_SPARSE_VECTOR_H

#include <memory>
#include <algorithm>
#include <limits>
#include <boost/operators.hpp>
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

(3) This sequence of index-value pairs is not kept sorted.
Most operations do however require that there is no duplicate index but since
this is too expensive to enforce in general, we leave it to user code
- either to  enforce it by construction,
- or to call the member function "sort_indices" when appropriate.
Any method requiring the abscence of duplicate will be annotated with:
precondition: no duplicate.

Implementation note:
The C++ standard rules that the private types and members of a class are not
accessible to its nested classes. A defect report (issue 45, [1]) has however
overturned that decision. It does not have the status "TC1" and is therefore not
part of the C++ standard yet. However, among the compilers supported by the cctbx,
VS 7.1 and 8.1, as well as all the GNU's have implemented the recommendation of
issue 45, i.e. that nested classes have access to all members of the outer
class. However cxx on Tru64 is true to the standard and does not implement it.
Therefore all inner classes have been made friends of the outer class vector.

[1] http://www.open-std.org/jtc1/sc22/wg21/docs/cwg_defects.html
*/
template<typename T>
class vector
{
  public:
    typedef T value_type;
    typedef std::size_t index_type;
    typedef std::ptrdiff_t index_difference_type;

  private:
    struct element {
      index_type index;
      value_type value;
      element(index_type i, value_type x) : index(i), value(x) {}

    };

    struct indexes_greater_than {
      bool operator()(const element& e, const element& f) {
        return e.index > f.index;
      }
    };
    friend struct indexes_greater_than;

    struct indexes_equal {
      bool operator()(const element& e, const element& f) {
        return e.index == f.index;
      }
    };
    friend struct indexes_equal;

    struct index_equal {
      index_type i;
      index_equal(index_type j) : i(j)
      {}
      bool operator()(const element& e) {
        return e.index == i;
      }
    };
    friend struct index_equal;

    typedef af::shared<element> container_type;

    container_type elements;
    index_type _size;

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
          return p->value;
        }

        bool operator==(const_iterator const& i) const {
          return p == i.p;
        }

        const_iterator& operator++() {
          p++;
          return *this;
        }

        index_type index() const {
          return p->index;
        }
    };
    friend class const_iterator;

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

        value_type& operator*() {
          return p->value;
        }

        bool operator==(iterator const& i) const {
          return p == i.p;
        }

        iterator& operator++() {
          p++;
          return *this;
        }

        index_type index() const {
          return p->index;
        }
    };
    friend class iterator;

    /// A reference to an element of given index
    /** This the type of object returned by v[i] for a vector v and an index i
    */
    class element_reference
    {
      protected:
        vector& v;
        index_type i;

      public:
        /// Construct the reference to the element of index j in u
        element_reference(vector& u, index_type j) : v(u), i(j)
        {}

        /* Without the destructor, the constructor would not be inlined
        by g++ 4.0 */
        ~element_reference()
        {}

        /// Triggered by using the value of v[i]
        operator value_type() {
          vector const &cv = v;
          typename container_type::const_reverse_iterator p = std::find_if(
            cv.elements.rbegin(), cv.elements.rend(), index_equal(i));
          return p != cv.elements.rend() ? p->value : 0;
        }

        /// Triggered by an assignment v[i] = ...
        value_type operator=(value_type x) {
          v.elements.push_back(element(i,x));
          return x;
        }
    };
    friend class element_reference;

    /// Construct a vector of size 0
    vector() {}

    /// Construct a zero vector of size n
    vector(index_type n) : _size(n) {}

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
      return _size;
    }

    /// Whether there is no potential non-zero elements
    bool is_structurally_zero() const {
      return elements.size() == 0;
    }

    /// Subscripting
    /** Assignment v[i] = ... may introduce a duplicate index, a problem
    solved for the getter v[i] by returning the last record with the
    desired index.
    */
    element_reference operator[](index_type i) {
      return element_reference(*this, i);
    }

    /// Fill the given dense vector.
    /** precondition: no duplicate */
    void fill_dense_vector(af::shared<T>& w) const {
      SCITBX_ASSERT(w.size() == size())
                   ( w.size() )( size() );
      for(const_iterator q =  begin(); q != end(); q++) {
        w[q.index()] = *q;
      }
    }

    /// The dense vector corresponding to this
    /** precondition: no duplicate */
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
      return result;
    }

    /// Fill the given dense vector with a permutation of this
    /** precondition: no duplicate */
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

    /// Remove duplicate indices and sort records by increasing indices.
    /** The record which was input last is kept in case of duplicates.
    Return this object, for convenient chaining of operations. */
    vector const& sort_indices() const {
      container_type &elts = const_cast<container_type&>(elements);
      std::stable_sort(elts.begin(), elts.end(), indexes_greater_than());
      container_type new_elements(af::reserve(elts.size()));
      std::unique_copy(elts.rbegin(), elts.rend(),
                       std::back_inserter(new_elements), indexes_equal());
      const_cast<vector*>(this)->elements = new_elements;
      return *this;
    }

    /// Permute the elements of this, in place
    /** Return this object, for convenient chaining of operations */
    template<class PermutationType>
    vector& permute(PermutationType const& permutation) {
      SCITBX_ASSERT(size() == permutation.size())
                   ( size() )( permutation.size() );
      for (typename container_type::iterator p=elements.begin();
           p != elements.end(); p++)
      {
        p->index = permutation[p->index];
      }
      return *this;
    }
};

/// Element-wise comparison of a and b with the given absolute tolerance
template<class T>
bool approx_equal(vector<T> const& a, vector<T> const& b,
                  T tol=std::numeric_limits<T>::epsilon())
{
    SCITBX_ASSERT(a.size() == b.size())
                ( a.size() )( b.size() );
    a.sort_indices();
    b.sort_indices();
    typename vector<T>::const_iterator p = a.begin(),
                                       q = b.begin();
    for(;;) {
      if (p == a.end() && q == b.end()) break;
      if (*p != 0 && *q != 0) {
        if (std::abs(*p - *q) > tol) return false;
        if (p != a.end()) p++;
        if (q != b.end()) q++;
      }
      else {
        if (*p == 0) p++;
        if (*q == 0) q++;
      }
    }
    return p == a.end() && q == b.end();
}


}}

#endif
