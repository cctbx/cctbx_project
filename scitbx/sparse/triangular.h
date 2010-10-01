#ifndef SCITBX_SPARSE_TRIANGULAR_H
#define SCITBX_SPARSE_TRIANGULAR_H

#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/accessors/packed_matrix.h>

namespace scitbx { namespace sparse {

  /// A proxy for the upper diagonal of a sparse matrix
  /** It is aimed at operating efficiently on a dense matrix.
   */
  template <class T>
  struct upper_diagonal
    : af::expression< upper_diagonal<T> >,
      operators::matrix_operating_on_dense_matrix<T, operators::upper_diagonal,
                                                  upper_diagonal<T> >
  {
    /// Minimal subset of matrix<T> interface
    //@{
    typedef typename matrix<T>::const_row_iterator const_row_iterator;
    typedef typename matrix<T>::column_type        column_type;
    typedef typename matrix<T>::index_type         index_type;
    typedef typename matrix<T>::value_type         value_type;

    int n_rows() const { return a.n_rows(); }
    int n_cols() const { return a.n_cols(); }
    column_type const &col(int j) const { return a.col(j); }
    column_type &col(int j) { return a.col(j); }
    void compact() const { a.compact(); }
    //@}

    matrix<T> const &a;
    upper_diagonal(matrix<T> const &a)
    : a(a)
    {
      SCITBX_ASSERT(n_rows() == n_cols())(n_rows())(n_cols());
    }

    af::packed_u_accessor
    expression_accessor(af::packed_u_accessor const &) const
    {
      return af::packed_u_accessor(n_rows());
    }

    /// Assign this to the given reference to a dense matrix
    /** This enables af::ref<T, af::packed_u_accessor> b = a
     for any sparse::matrix<T> a
     */
    void assign_to(af::ref<T, af::packed_u_accessor> const &b) const {
      this->operate_on(operators::assign<T>(), b);
    }

    /// Add this to the given reference to a dense matrix
    /** This enables af::ref<T, af::packed_u_accessor> b += a
     for any sparse::matrix<T> a
     */
    void add_to(af::ref<T, af::packed_u_accessor> const &b) const {
      this->operate_on(operators::plus_equal<T>(), b);
    }

    /// Substract from the given reference to a dense matrix
    /** This enables af::ref<T, af::packed_u_accessor> b -= a
     for any sparse::matrix<T> a
     */
    void substract_from(af::ref<T, af::packed_u_accessor> const &b) const {
      this->operate_on(operators::minus_equal<T>(), b);
    }

  };

  /// The usual convenience factor function deducing the template argument
  template <class T>
  upper_diagonal<T> upper_diagonal_of(matrix<T> const &a) {
    return upper_diagonal<T>(a);
  }

}}

#endif // GUARD
