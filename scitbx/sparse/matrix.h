#ifndef SCITBX_SPARSE_MATRIX_H
#define SCITBX_SPARSE_MATRIX_H

#include <algorithm>
#include <functional>
#include <vector>
#include <scitbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/sparse/vector.h>

namespace scitbx { namespace sparse {

/// A sparse matrix, represented by a sequence of sparse columns
/** All linear operations are therefore performed using column version of the
 relevant algorithms, taking great care of never touching structurally zero
 elements.

 Assignments a(i, j) = ..., a(i,j) += ..., and a(i,j) -= ... are delegated
 to the corresponding assignments c[i] = ..., c[i] += ... and c[i] -= ...
 where c[i] is the j-th column col(j), an instance of sparse::vector, to
 the documentation of which the reader is referred to.
 It should be noticed that no range checking on i is performed for those
 assignments, but that fetching a(i,j) with i >= n_rows() will always return 0.
*/
template<class T>
class matrix
{
private:
  /* Let's take advantage of scitbx reference counted array semantic...
  C.f., for example, member function "transpose"
  */
  typedef af::shared< vector<T> > container_type;
  typedef af::ref< vector<T> > ref_type;

public:
  typedef T value_type;
  typedef vector<T> column_type;
  typedef typename vector<T>::index_type row_index;
  typedef typename container_type::size_type column_index;
  typedef typename vector<T>::iterator row_iterator;
  typedef typename vector<T>::const_iterator const_row_iterator;

public:
  /// Construct a zero matrix with the given number of rows and columns
  matrix(row_index rows, column_index cols)
    : n_rows_(rows), columns(af::reserve(rows))
  {
    for (column_index j=0; j < cols; j++) {
      columns.push_back(column_type(rows));
    }
    column = columns.ref();
  }

  /// The i-th column
  vector<T>& col(column_index i) {
    return column[i];
  }

  /// The i-th column
  vector<T> const& col(column_index i) const {
    return column[i];
  }

  /// Subscripting
  /** This pulls out the column j and then delegates subscripting of row
  index i to class vector, which the readers is referred to.
  */
  typename vector<T>::element_reference
  operator()(row_index i, column_index j) {
    return column[j][i];
  }

  /// Subscripting
  /** This pulls out the column j and then delegates subscripting of row
    index i to class vector, which the readers is referred to.
    */
  const typename vector<T>::element_const_reference
  operator()(row_index i, column_index j) const {
    return column[j][i];
  }

  /// Whether the element (i,j) is a structural zero
  bool is_structural_zero(row_index i, column_index j) {
    return column[j].is_structural_zero(i);
  }

  /// Number of columns
  column_index n_cols() const {
    return column.size();
  }

  /// Number of rows
  row_index n_rows() const {
    return n_rows_;
  }

  /// Whether all elements below the diagonal are zero
  bool is_upper_triangular() const {
    for (column_index j=0; j < n_cols(); j++) {
      for (const_row_iterator p = col(j).begin(); p != col(j).end(); p++) {
        if (p.index() > j && *p != 0) return false;
      }
    }
    return true;
  }

  /// Whether all elements above the diagonal are zero and the diagonal is 1
  bool is_unit_lower_triangular() const {
    for (column_index j=0; j < n_cols(); j++) {
      for (const_row_iterator p = col(j).begin(); p != col(j).end(); p++) {
        if (p.index() < j && *p != 0) return false;
        else if (p.index() == j && *p != 1) return false;
      }
    }
    return true;
  }

  /// Number of non-zero elements
  index_type non_zeroes() const {
    index_type result = 0;
    for (int j=0; j<n_cols(); ++j) result += col(j).non_zeroes();
    return result;
  }

  /// A copy of this matrix, copying elements
  matrix deep_copy() const {
    matrix result(n_rows(), n_cols());
    for (column_index j=0; j < n_cols(); j++) {
      result.column[j] = column[j].deep_copy();
    }
    return result;
  }

  /// Transpose of this
  matrix transpose() const {
    matrix result(n_cols(), n_rows());
    for (column_index j=0; j < n_cols(); j++) {
      for(const_row_iterator p = col(j).begin(); p != col(j).end(); p++) {
        result(j, p.index()) = *p;
      }
    }
    return result;
  }

  /// Sort and remove duplicate indices in all column
  /** C.f. the member function of same name in class scitbx::sparse::vector
   In particular, any element which may have been introduced into the matrix
   that has a row index i >= n_rows() is pruned.
   */
  void compact() const {
    for (column_index j=0; j < n_cols(); j++) col(j).compact();
  }

  /// Permute the rows, in place
  template<class PermutationType>
  matrix& permute_rows(PermutationType const& permutation) {
    SCITBX_ASSERT(n_rows() == permutation.size())
                 ( n_rows() )( permutation.size() );
    for (column_index j=0; j < n_cols(); j++) col(j).permute(permutation);
    return *this;
  }

  /// This times sparse column vector
  vector<T> operator*(vector<T> const& v) const {
    SCITBX_ASSERT(n_cols() == v.size())
                 ( n_cols() )( v.size() );
    vector<T> w(n_rows());
    for (const_row_iterator q=v.begin(); q!=v.end(); ++q) {
      column_index j = q.index();
      value_type v_j = *q;
      for (const_row_iterator p=col(j).begin(); p != col(j).end(); ++p) {
        row_index i = p.index();
        value_type a_ij = *p;
        w[i] += a_ij * v_j;
      }
    }
    w.compact();
    return w;
  }

  typedef typename vector<T>::dense_vector_const_ref dense_vector_const_ref;
  typedef af::shared<T> dense_vector;

  /// This times dense column vector
  /** Our returning a dense vector here is most appropriate when
   there are few zero rows.
   */
  dense_vector operator*(dense_vector_const_ref const &v) const {
    SCITBX_ASSERT(n_cols() == v.size())
    ( n_cols() )( v.size() );
    dense_vector w(n_rows());
    for (int j=0; j < n_cols(); ++j) {
      for (const_row_iterator p=col(j).begin(); p != col(j).end(); ++p) {
        row_index i = p.index();
        value_type a_ij = *p;
        w[i] += a_ij * v[j];
      }
    }
    return w;
  }

  /// Dense row vector times matrix
  /** Our returning a dense vector here is most appropriate when
   there are few zero columns.
   */
  friend
  dense_vector operator*(dense_vector_const_ref const &u, matrix const &a) {
    dense_vector result(a.n_cols(), af::init_functor_null<value_type>());
    for(column_index j=0; j < a.n_cols(); ++j) result[j] = u*a.col(j);
    return result;
  }

  /// Matrix times matrix
  friend
  matrix operator*(matrix const& a, matrix const& b) {
    SCITBX_ASSERT(a.n_cols() == b.n_rows())
                 ( a.n_cols() )( b.n_rows() );
    matrix c(a.n_rows(), b.n_cols());
    for (column_index j=0; j < c.n_cols(); j++) c.col(j) = a * b.col(j);
    return c;
  }

  typedef af::versa<value_type, af::packed_u_accessor>
          symmetric_matrix_t;
  typedef af::const_ref<value_type, af::packed_u_accessor>
          symmetric_matrix_const_ref_t;
  typedef af::ref<value_type, af::packed_u_accessor>
          symmetric_matrix_ref_t;

  /// B^T A B where B is this matrix and A is a symmetric dense matrix
      /** This is useful for the interplay between change of variable and
      least-squares covariance matrix. If x and y are parameter vectors
      related by of respective size p < n:
      \f[
          \nabla_x = B \nabla_y,
      \f]
      where B is a p x n matrix, or equivalently
      \f[
          \delta y = B^T \delta x,
      \f]
       then the covariance matrices \f$V_y\f$ and \f$V_x\f$ for respectively
       x and y are related by
      \f[
          V_y = B^T V_x B. (1)
      \f]
      In practice,
      \f[
          B = \left[ \frac{\partial y_j}{\partial x_i} \right]_{ij}
      \f]
      has sparse columns because \f$y_j\f$ depends only on a few \f$x_i\f$,
      which makes the storage scheme adopted here a perfect match memory-wise.
      Furthermore, the computation of (1) turns out to be well adapted
      to that storage scheme too.
   */
  symmetric_matrix_t
  this_transpose_times_symmetric_times_this(
    symmetric_matrix_const_ref_t const &a) const
  {
    compact();

    // This verbosity seems necessary to please gcc 3.2
    af::packed_u_accessor sym_n_x_n(n_cols());
    af::init_functor_null<value_type> dont_init;
    symmetric_matrix_t result(sym_n_x_n, dont_init);
    symmetric_matrix_ref_t c = result.ref();

    int n = c.accessor().n;
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
      c(i,j) = 0;
      for (const_row_iterator p = col(i).begin(); p != col(i).end(); ++p) {
        row_index k = p.index();
        value_type b_ki = *p;
        value_type s = 0;
        for (const_row_iterator q = col(j).begin(); q != col(j).end(); ++q) {
          row_index l = q.index();
          value_type b_lj = *q;
          if (k <= l) s += a(k,l)*b_lj;
          else        s += a(l,k)*b_lj;
        }
        c(i,j) += b_ki*s;
      }
    }
    return result;
  }

private:
  typedef typename std::vector<row_index>::const_iterator const_row_idx_iter;
  container_type columns;
  ref_type column;
  row_index n_rows_;
};

}} // namespace scitbx::sparse


#endif
