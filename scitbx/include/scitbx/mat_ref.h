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

  typedef
    af::grid<2, af::c_index_1d_calculator<2>, af::tiny<std::size_t, 2> >
      mat_grid;

  template <typename NumType, typename AccessorType = mat_grid>
  class mat_ref : public af::ref<NumType, AccessorType>
  {
    public:
      typedef AccessorType accessor_type;
      typedef typename af::ref<NumType, AccessorType> base_type;
      typedef typename accessor_type::index_type::value_type size_type;

      mat_ref() {}

      mat_ref(NumType* begin, accessor_type const& grid)
      : base_type(begin, grid)
      {}

      mat_ref(NumType* begin, size_type n_rows, size_type n_columns)
      : base_type(begin, accessor_type(n_rows, n_columns))
      {}

      accessor_type grid() const { return this->accessor(); }

      size_type n_rows() const { return this->accessor()[0]; }

      size_type n_columns() const { return this->accessor()[1]; }

      bool is_same_grid(mat_ref const& other) const
      {
        if (n_rows() != other.n_rows()) return false;
        if (n_columns() != other.n_columns()) return false;
        return true;
      }

      //! Access elements with 2-dimensional indices.
      NumType const&
      operator()(size_type r, size_type c) const
      {
        return this->begin()[r * n_columns() + c];
      }

      //! Access elements with 2-dimensional indices.
      NumType&
      operator()(size_type r, size_type c)
      {
        return this->begin()[r * n_columns() + c];
      }

      //! Test for square matrix.
      bool is_square() const { return n_rows() == n_columns(); }

      //! Test for diagonal matrix.
      bool is_diagonal() const;

      //! Swap two rows in place.
      void
      swap_rows(size_type i1, size_type i2)
      {
        std::swap_ranges(&(*this)(i1,0), &(*this)(i1+1,0), &(*this)(i2,0));
      }

      //! Swap two columns in place.
      void
      swap_columns(size_type i1, size_type i2)
      {
        for(size_type ir=0;ir<n_rows();ir++) {
          std::swap((*this)(ir,i1), (*this)(ir,i2));
        }
      }

      //! Sets diagonal matrix. Requires a square matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_diagonal(NumType const& d);

      //! Sets identity matrix. Requires a square matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_identity()
      {
        set_diagonal(1);
      }

      //! Transposes matrix in place.
      void transpose_in_place();
  };

  // non-inline member function
  template <typename NumType, typename AccessorType>
  bool
  mat_ref<NumType, AccessorType>
  ::is_diagonal() const
  {
    if (!is_square()) return false;
    for (size_type ir=0;ir<n_rows();ir++)
      for (size_type ic=0;ic<n_columns();ic++)
        if (ir != ic && (*this)(ir,ic)) return false;
    return true;
  }

  // non-inline member function
  template <typename NumType, typename AccessorType>
  void
  mat_ref<NumType, AccessorType>
  ::set_diagonal(NumType const& d)
  {
   SCITBX_ASSERT(is_square());
   this->fill(0);
   for(size_type i=0;i<n_rows();i++) (*this)(i,i) = d;
  }

  // non-inline member function
  template <typename NumType, typename AccessorType>
  void
  mat_ref<NumType, AccessorType>
  ::transpose_in_place()
  {
    if (is_square()) {
      for (size_type ir=0;ir<n_rows();ir++)
        for (size_type ic=ir+1;ic<n_columns();ic++)
          std::swap((*this)(ir, ic), (*this)(ic, ir));
    }
    else {
      std::vector<NumType> mt_buffer(this->size());
      mat_ref mt(&*mt_buffer.begin(), n_columns(), n_rows());
      for (size_type ir=0;ir<n_rows();ir++)
        for (size_type ic=0;ic<n_columns();ic++)
          mt(ic, ir) = (*this)(ir, ic);
      std::copy(mt.begin(), mt.end(), this->begin());
      this->m_accessor = mt.accessor();
    }
  }

  //! Test equality.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    mat_ref<NumType, AccessorType> const& lhs,
    mat_ref<NumType, AccessorType> const& rhs)
  {
    if (!lhs.is_same_grid(rhs)) return false;
    return !af::order(lhs, rhs);
  }

  //! Test equality. True if all elements of lhs == rhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    mat_ref<NumType, AccessorType> const& lhs,
    NumType const& rhs)
  {
    return !af::order(lhs, rhs);
  }

  //! Test equality. True if all elements of rhs == lhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator==(
    NumType const& lhs,
    mat_ref<NumType, AccessorType> const& rhs)
  {
    return !af::order(lhs, rhs);
  }

  //! Test inequality.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator!=(
    mat_ref<NumType, AccessorType> const& lhs,
    mat_ref<NumType, AccessorType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Test inequality. True if any element of lhs != rhs.
  template <typename NumType, typename AccessorType>
  inline
  bool
  operator!=(
    mat_ref<NumType, AccessorType> const& lhs,
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
    mat_ref<NumType, AccessorType> const& rhs)
  {
    return !(lhs == rhs);
  }

  //! Element-wise in-place addition.
  template <typename NumType, typename AccessorType>
  inline
  mat_ref<NumType, AccessorType>&
  operator+=(
    mat_ref<NumType, AccessorType>& lhs,
    mat_ref<NumType, AccessorType> const& rhs)
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
  mat_ref<NumType, AccessorType>&
  operator+=(
    mat_ref<NumType, AccessorType>& lhs,
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
  mat_ref<NumType, AccessorType>&
  operator-=(
    mat_ref<NumType, AccessorType>& lhs,
    mat_ref<NumType, AccessorType> const& rhs)
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
  mat_ref<NumType, AccessorType>&
  operator-=(
    mat_ref<NumType, AccessorType>& lhs,
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
  mat_ref<NumType, AccessorType>&
  operator*=(
    mat_ref<NumType, AccessorType>& lhs,
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
  mat_ref<NumType, AccessorType>&
  operator/=(
    mat_ref<NumType, AccessorType>& lhs,
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
  mat_ref<NumType, AccessorType>&
  operator%=(
    mat_ref<NumType, AccessorType>& lhs,
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
    mat_ref<NumTypeA, AccessorTypeA> const& a,
    mat_ref<NumTypeB, AccessorTypeB> const& b,
    mat_ref<NumTypeAB, AccessorTypeAB>& ab)
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
