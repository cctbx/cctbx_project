/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_MAT_REF_H
#define SCITBX_MAT_REF_H

#include <scitbx/error.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/matrix_multiply.h>
#include <vector>

namespace scitbx {

  class mat_grid : public af::tiny<std::size_t, 2>
  {
    public:
      typedef af::tiny<std::size_t, 2> index_type;
      typedef index_type::value_type index_value_type;

      mat_grid() : index_type(0,0) {}

      mat_grid(index_type const& n) : index_type(n) {}

      mat_grid(index_value_type const& n0, index_value_type const& n1)
      : index_type(n0, n1)
      {}

      std::size_t size_1d() const { return elems[0] * elems[1]; }

      std::size_t
      operator()(index_value_type const& r, index_value_type const& c) const
      {
        return r * elems[1] + c;
      }
  };

  template <typename NumType, typename AccessorType = mat_grid>
  class mat_const_ref : public af::const_ref<NumType, AccessorType>
  {
    public:
      typedef AccessorType accessor_type;
      typedef typename af::const_ref<NumType, AccessorType> base_type;
      typedef typename accessor_type::index_value_type index_value_type;

      mat_const_ref() {}

      mat_const_ref(const NumType* begin, accessor_type const& grid)
      : base_type(begin, grid)
      {}

      mat_const_ref(const NumType* begin, index_value_type const& n_rows,
                                          index_value_type const& n_columns)
      : base_type(begin, accessor_type(n_rows, n_columns))
      {}

      accessor_type
      grid() const { return this->accessor(); }

      index_value_type const&
      n_rows() const { return this->accessor()[0]; }

      index_value_type const&
      n_columns() const { return this->accessor()[1]; }

      //! Tests for square matrix.
      bool
      is_square() const { return n_rows() == n_columns(); }

      bool
      is_same_grid(mat_const_ref const& other) const
      {
        if (n_rows() != other.n_rows()) return false;
        if (n_columns() != other.n_columns()) return false;
        return true;
      }

      //! Accesses elements with 2-dimensional indices.
      NumType const&
      operator()(index_value_type const& r, index_value_type const& c) const
      {
        return this->begin()[this->accessor()(r, c)];
      }

      //! Tests for diagonal matrix.
      bool
      is_diagonal() const;
  };

  // non-inline member function
  template <typename NumType, typename AccessorType>
  bool
  mat_const_ref<NumType, AccessorType>
  ::is_diagonal() const
  {
    if (!is_square()) return false;
    for (index_value_type ir=0;ir<n_rows();ir++)
      for (index_value_type ic=0;ic<n_columns();ic++)
        if (ir != ic && (*this)(ir,ic)) return false;
    return true;
  }

  template <typename NumType, typename AccessorType = mat_grid>
  class mat_ref : public mat_const_ref<NumType, AccessorType>
  {
    public:
      typedef AccessorType accessor_type;
      typedef mat_const_ref<NumType, AccessorType> base_type;
      typedef typename accessor_type::index_value_type index_value_type;

      mat_ref() {}

      mat_ref(NumType* begin, accessor_type const& grid)
      : base_type(begin, grid)
      {}

      mat_ref(NumType* begin, index_value_type n_rows,
                              index_value_type n_columns)
      : base_type(begin, accessor_type(n_rows, n_columns))
      {}

      NumType*
      begin() const { return const_cast<NumType*>(this->begin_); }

      NumType*
      end() const { return const_cast<NumType*>(this->end_); }

      NumType&
      front() const { return begin()[0]; }

      NumType&
      back() const { return end()[-1]; }

      NumType&
      operator[](index_value_type const& i) const { return begin()[i]; }

      NumType&
      at(index_value_type const& i) const
      {
        if (i >= this->size()) af::throw_range_error();
        return begin()[i];
      }

      mat_ref const&
      fill(NumType const& x) const
      {
        std::fill(begin(), end(), x);
        return *this;
      }

      //! Accesses elements with 2-dimensional indices.
      NumType&
      operator()(index_value_type const& r, index_value_type const& c) const
      {
        return this->begin()[this->accessor()(r, c)];
      }

      //! Swaps two rows in place.
      void
      swap_rows(index_value_type const& i1, index_value_type const& i2) const
      {
        std::swap_ranges(&(*this)(i1,0), &(*this)(i1+1,0), &(*this)(i2,0));
      }

      //! Swaps two columns in place.
      void
      swap_columns(index_value_type const& i1,
                   index_value_type const& i2) const
      {
        for(index_value_type ir=0;ir<this->n_rows();ir++) {
          std::swap((*this)(ir,i1), (*this)(ir,i2));
        }
      }

      //! Sets diagonal matrix. Requires a square matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_diagonal(NumType const& d) const;

      //! Sets identity matrix. Requires a square matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_identity() const
      {
        set_diagonal(1);
      }

      void transpose_square_in_place() const
      {
        SCITBX_ASSERT(this->is_square());
        for (index_value_type ir=0;ir<this->n_rows();ir++)
          for (index_value_type ic=ir+1;ic<this->n_columns();ic++)
            std::swap((*this)(ir, ic), (*this)(ic, ir));
      }

      //! Transposes matrix in place.
      void transpose_in_place();
  };

  // non-inline member function
  template <typename NumType, typename AccessorType>
  void
  mat_ref<NumType, AccessorType>
  ::set_diagonal(NumType const& d) const
  {
   SCITBX_ASSERT(this->is_square());
   this->fill(0);
   for(index_value_type i=0;i<this->n_rows();i++) (*this)(i,i) = d;
  }

  // non-inline member function
  template <typename NumType, typename AccessorType>
  void
  mat_ref<NumType, AccessorType>
  ::transpose_in_place()
  {
    if (this->is_square()) {
      for (index_value_type ir=0;ir<this->n_rows();ir++)
        for (index_value_type ic=ir+1;ic<this->n_columns();ic++)
          std::swap((*this)(ir, ic), (*this)(ic, ir));
    }
    else {
      std::vector<NumType> mt_buffer(this->size());
      mat_ref mt(&*mt_buffer.begin(), this->n_columns(), this->n_rows());
      for (index_value_type ir=0;ir<this->n_rows();ir++)
        for (index_value_type ic=0;ic<this->n_columns();ic++)
          mt(ic, ir) = (*this)(ir, ic);
      std::copy(mt.begin(), mt.end(), this->begin());
      this->accessor_ = mt.accessor();
      this->init();
    }
  }

  //! Tests equality.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    mat_const_ref<NumType, AccessorType> const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    if (!lhs.is_same_grid(rhs)) return false;
    return lhs.all_eq(rhs);
  }

  //! Tests equality. True if all elements of lhs == rhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    mat_const_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    return lhs.all_eq(rhs);
  }

  //! Test equality. True if all elements of rhs == lhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    NumType const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    return rhs.all_eq(lhs);
  }

  //! Test inequality.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator!=(
    mat_const_ref<NumType, AccessorType> const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator!=(
    mat_const_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of rhs != lhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator!=(
    NumType const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise in-place addition.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator+=(
    mat_ref<NumType, AccessorType> const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    SCITBX_ASSERT(lhs.is_same_grid(rhs));
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] += rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place addition.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator+=(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] += rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator-=(
    mat_ref<NumType, AccessorType> const& lhs,
    mat_const_ref<NumType, AccessorType> const& rhs)
  {
    SCITBX_ASSERT(lhs.is_same_grid(rhs));
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] -= rhs[i];
    }
    return lhs;
  }

  //! Element-wise in-place difference.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator-=(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] -= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place multiplication.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator*=(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] *= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place division.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator/=(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] /= rhs   ;
    }
    return lhs;
  }

  //! Element-wise in-place modulus operation.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType> const&
  operator%=(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    for(std::size_t i=0;i<lhs.size();i++) {
      lhs[i] %= rhs   ;
    }
    return lhs;
  }

  //! Matrix multiplication.
  template <typename NumTypeA, typename AccessorTypeA,
            typename NumTypeB, typename AccessorTypeB,
            typename NumTypeAB, typename AccessorTypeAB>
  inline
  void
  multiply(
    mat_const_ref<NumTypeA, AccessorTypeA> const& a,
    mat_const_ref<NumTypeB, AccessorTypeB> const& b,
    mat_ref<NumTypeAB, AccessorTypeAB> const& ab)
  {
    SCITBX_ASSERT(a.n_columns() == b.n_rows());
    SCITBX_ASSERT(ab.n_rows() == a.n_rows());
    SCITBX_ASSERT(ab.n_columns() == b.n_columns());
    matrix_multiply(a.begin(), b.begin(),
                    a.n_rows(), a.n_columns(), b.n_columns(),
                    ab.begin());
  }

} // namespace scitbx

#endif // SCITBX_MAT_REF_H
