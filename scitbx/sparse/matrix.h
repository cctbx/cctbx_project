#ifndef SCITBX_SPARSE_MATRIX_H
#define SCITBX_SPARSE_MATRIX_H

#include <algorithm>
#include <functional>
#include <vector>
#include <scitbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/util.h>
#include <scitbx/sparse/operators.h>

namespace scitbx { namespace sparse {

template <class T>
class matrix;

template <class T>
struct matrix_times_dense_vector
  : af::expression< matrix_times_dense_vector<T> >
{
  typedef T value_type;

  matrix_times_dense_vector(matrix<T> const &a, af::const_ref<T> const &v)
  : a(a), v(v)
  {}

  std::size_t size() const { return a.n_rows(); }

  void assign_to(af::ref<T> const &w) const;

  matrix<T> const &a;
  af::const_ref<T> const &v;
};

template <class T>
struct matrix_transpose_times_dense_vector
  : af::expression< matrix_transpose_times_dense_vector<T> >
{
  typedef T value_type;

  matrix_transpose_times_dense_vector(matrix<T> const &a,
                                      af::const_ref<T> const &v)
  : a(a), v(v)
  {}

  std::size_t size() const { return a.n_cols(); }

  void assign_to(af::ref<T> const &w) const;

  matrix<T> const &a;
  af::const_ref<T> const &v;
};



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

 The columns have a copy semantic.
*/
template<class T>
class matrix
  : public boost::equality_comparable< matrix<T> >,
    public af::expression< matrix<T> >,
    public operators::matrix_operating_on_dense_matrix<T, operators::whole,
                                                       matrix<T> >
{
public:
  typedef T value_type;

  /* Use std::vector as the underlying container for the sparse columns
     in an effort to reduce the memory overhead.
     C.f. flex.show_sizes_double for the rationales behind not using af::shared
   */
  typedef vector<T, copy_semantic_vector_container> column_type;
  typedef typename column_type::index_type index_type;
  typedef typename column_type::iterator row_iterator;
  typedef typename column_type::const_iterator const_row_iterator;

public:
  /// Construct uninitialised object
  /** The only safe operation is to assign a well initialised sparse::matrix
   */
  matrix() {}

  /// Construct a zero matrix with the given number of rows and columns
  matrix(index_type rows, index_type cols)
    : n_rows_(rows), columns(af::reserve(cols))
  {
    for (index_type j=0; j < cols; j++) {
      columns.push_back(column_type(rows));
    }
    column = columns.ref();
  }

  /// The i-th column
  column_type& col(index_type i) {
    return column[i];
  }

  /// The i-th column
  column_type const& col(index_type i) const {
    return column[i];
  }

  /// Subscripting
  /** This pulls out the column j and then delegates subscripting of row
  index i to class vector, which the readers is referred to.
  */
  typename column_type::element_reference
  operator()(index_type i, index_type j) {
    return column[j][i];
  }

  /// Subscripting
  /** This pulls out the column j and then delegates subscripting of row
    index i to class vector, which the readers is referred to.
    */
  const typename column_type::element_const_reference
  operator()(index_type i, index_type j) const {
    return column[j][i];
  }

  /// Assign the given matrix to the block starting at row i and column j
  void assign_block(matrix const &b, int i, int j) {
    SCITBX_ASSERT(i + b.n_rows() <= n_rows())(i)(b.n_rows())(n_rows());
    SCITBX_ASSERT(j + b.n_cols() <= n_cols())(j)(b.n_cols())(n_cols());
    for (int jj=0; jj<b.n_cols(); ++jj) {
      for (const_row_iterator p=b.col(jj).begin(); p != b.col(jj).end(); ++p)
      {
        index_type ii = p.index();
        (*this)(i + ii, j + jj) = *p;
      }
    }
    compact();
  }

  /// Assign the given matrix to the block starting at row i and column j
  /** Zero elements in b results in structural zeroes in the assign sparse block
   */
  void assign_block(af::const_ref<T, af::mat_grid> const &b, int i, int j) {
    SCITBX_ASSERT(i + b.n_rows() <= n_rows())(i)(b.n_rows())(n_rows());
    SCITBX_ASSERT(j + b.n_columns() <= n_cols())(j)(b.n_columns())(n_cols());
    for (int jj=0; jj<b.n_columns(); ++jj) for (int ii=0; ii<b.n_rows(); ++ii)
    {
      if (b(ii, jj) != T(0)) (*this)(i + ii, j + jj) = b(ii, jj);
    }
    compact();
  }

  /// Whether the element (i,j) is a structural zero
  bool is_structural_zero(index_type i, index_type j) {
    return column[j].is_structural_zero(i);
  }

  /// Number of columns
  index_type n_cols() const {
    return column.size();
  }

  /// Number of rows
  index_type n_rows() const {
    return n_rows_;
  }

  /// Whether all elements below the diagonal are zero
  bool is_upper_triangular() const {
    for (index_type j=0; j < n_cols(); j++) {
      for (const_row_iterator p = col(j).begin(); p != col(j).end(); p++) {
        if (p.index() > j && *p != 0) return false;
      }
    }
    return true;
  }

  /// Whether all elements above the diagonal are zero and the diagonal is 1
  bool is_unit_lower_triangular() const {
    for (index_type j=0; j < n_cols(); j++) {
      for (const_row_iterator p = col(j).begin(); p != col(j).end(); p++) {
        if (p.index() < j && *p != 0) return false;
        else if (p.index() == j && *p != 1) return false;
      }
    }
    return true;
  }

  /// Equality
  bool operator==(matrix const &other) const {
    return std::equal(columns.begin(), columns.end(), other.columns.begin());
  }

  /// Number of non-zero elements
  index_type non_zeroes() const {
    index_type result = 0;
    for (int j=0; j<n_cols(); ++j) result += col(j).non_zeroes();
    return result;
  }

  /// A copy of this matrix, copying elements
  matrix clone() const {
    matrix result(n_rows(), n_cols());
    for (index_type j=0; j < n_cols(); j++) {
      result.column[j] = column[j];
    }
    return result;
  }

  /// The matrix whose columns are the given (sub)set of this matrix columns
  /** If the elements of p are a permutation of the range 0..n-1 where n is
      the number of columns, then this is a column-permuted variant of this
      matrix.
   */
  matrix select_columns(af::const_ref<index_type> const &p) const {
    SCITBX_ASSERT(p.size() <= n_cols())(p.size())(n_cols());
    matrix result(n_rows(), p.size());
    for (index_type j=0; j<p.size(); ++j) result.col(j) = this->col(p[j]);
    return result;
  }

  /// Transpose of this
  matrix transpose() const {
    matrix result(n_cols(), n_rows());
    for (index_type j=0; j < n_cols(); j++) {
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
    for (index_type j=0; j < n_cols(); j++) col(j).compact();
  }

  /// Specify whether the matrix shall be considered compacted or not
  /** set_compact(true) is very dangerous: it shall only be used
      along with algorithms provingly building matrices with increasing
      unique row indices.
   */
  void set_compact(bool f) const {
    for (index_type j=0; j < n_cols(); j++) col(j).set_compact(f);
  }

  /// Permute the rows, in place
  template<class PermutationType>
  matrix& permute_rows(PermutationType const& permutation) {
    SCITBX_ASSERT(n_rows() == permutation.size())
                 ( n_rows() )( permutation.size() );
    for (index_type j=0; j < n_cols(); j++) col(j).permute(permutation);
    return *this;
  }

  /// This times sparse column vector
  template <class U, template<class> class C>
  vector<U,C> operator*(vector<U,C> const& v) const {
    SCITBX_ASSERT(n_cols() == v.size())
                 ( n_cols() )( v.size() );
    vector<U,C> w(n_rows());
    for (typename vector<U,C>::const_iterator q=v.begin(); q!=v.end(); ++q) {
      index_type j = q.index();
      value_type v_j = *q;
      for (const_row_iterator p=col(j).begin(); p != col(j).end(); ++p) {
        index_type i = p.index();
        value_type a_ij = *p;
        w[i] += a_ij * v_j;
      }
    }
    w.compact();
    return w;
  }

  typedef typename column_type::dense_vector_const_ref dense_vector_const_ref;
  typedef af::shared<T> dense_vector;

  /// This times dense column vector
  /** Our returning a dense vector here is most appropriate when
   there are few zero rows.
   */
  matrix_times_dense_vector<T> operator*(dense_vector_const_ref const &v) const {
    SCITBX_ASSERT(n_cols() == v.size())
    ( n_cols() )( v.size() );
    return matrix_times_dense_vector<T>(*this, v);
  }

  /// A^T v where A is this matrix
  matrix_transpose_times_dense_vector<T>
  transpose_times(dense_vector_const_ref const &v) const {
    SCITBX_ASSERT(n_rows() == v.size())( n_rows() )( v.size() );
    return matrix_transpose_times_dense_vector<T>(*this, v);
  }

  /// Dense row vector times matrix
  /** Our returning a dense vector here is most appropriate when
   there are few zero columns.
   */
  friend
  dense_vector operator*(dense_vector_const_ref const &u, matrix const &a) {
    dense_vector result(a.n_cols(), af::init_functor_null<value_type>());
    for(index_type j=0; j < a.n_cols(); ++j) result[j] = u*a.col(j);
    return result;
  }

  /// Accessor for a dense matrix that would correspond to this
  af::mat_grid expression_accessor(af::mat_grid const &proto) const {
    return af::mat_grid(n_rows(), n_cols());
  }

  /// Assign this to the given reference to a dense matrix
  /** This enables af::ref<T, af::mat_grid> b = a
      for any sparse::matrix<T> a
   */
  void assign_to(af::ref<T, af::mat_grid> const &b) const {
    this->operate_on(operators::assign<T>(), b);
  }

  /// Add this to the given reference to a dense matrix
  /** This enables af::ref<T, af::mat_grid> b += a
      for any sparse::matrix<T> a
   */
  void add_to(af::ref<T, af::mat_grid> const &b) const {
    this->operate_on(operators::plus_equal<T>(), b);
  }

  /// Substract from the given reference to a dense matrix
  /** This enables af::ref<T, af::mat_grid> b -= a
      for any sparse::matrix<T> a
   */
  void substract_from(af::ref<T, af::mat_grid> const &b) const {
    this->operate_on(operators::minus_equal<T>(), b);
  }

  /// Matrix times matrix
  friend
  matrix operator*(matrix const& a, matrix const& b) {
    SCITBX_ASSERT(a.n_cols() == b.n_rows())
                 ( a.n_cols() )( b.n_rows() );
    matrix c(a.n_rows(), b.n_cols());
    for (index_type j=0; j < c.n_cols(); j++) c.col(j) = a * b.col(j);
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
      of respective size p < n related by :
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
      This member function is well-adpated to the case where
      \f[
          B = \left[ \frac{\partial y_j}{\partial x_i} \right]_{ij}
      \f]
      has sparse columns, i.e. when each \f$y_j\f$ depends only
      on a few \f$x_i\f$.
   */
  symmetric_matrix_t
  this_transpose_times_symmetric_times_this(
    symmetric_matrix_const_ref_t const &a) const
  {
    SCITBX_ASSERT(a.accessor().n == n_rows());

    compact();

    af::packed_u_accessor sym_n_x_n(n_cols());
    symmetric_matrix_t result(sym_n_x_n);
    value_type *c = result.begin();
    int n = result.accessor().n;
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
      value_type &c_ij = *c++;
      for (const_row_iterator p = col(i).begin(); p != col(i).end(); ++p) {
        index_type k = p.index();
        value_type b_ki = *p;
        value_type s = 0;
        for (const_row_iterator q = col(j).begin(); q != col(j).end(); ++q) {
          index_type l = q.index();
          value_type b_lj = *q,
                     a_kl = k <= l ? a(k,l) : a(l,k);
          s += a_kl*b_lj;
        }
        c_ij += b_ki*s;
      }
    }
    return result;
  }

  /// B A B^T where B is this matrix and A is a dense symmetric matrix
  /** This has the same use as this_transpose_times_symmetric_times_this()
   except that this time
   \f[
      B = \left[ \frac{\partial y_i}{\partial x_j} \right]_{ij}
   \f]
   and therefore that this scheme is adpated to the case where for each
   variable \f$x_j\f$, there are only a few variable \f$y_i\f$ depending on it.
   */
  symmetric_matrix_t
  this_times_symmetric_times_this_transpose(
    symmetric_matrix_const_ref_t const &a) const
  {
    SCITBX_ASSERT(a.accessor().n == n_cols());

    compact();

    af::packed_u_accessor sym_m_x_m(n_rows());
    symmetric_matrix_t result(sym_m_x_m);
    symmetric_matrix_ref_t c = result.ref();
    int n = a.accessor().n;
    value_type const *a_ = a.begin();
    for (int k=0; k<n; ++k) {
      value_type a_kk = *a_++;
      for (const_row_iterator p = col(k).begin(); p != col(k).end(); ++p) {
        index_type i = p.index();
        value_type b_ik = *p;
        for (const_row_iterator q = p; q != col(k).end(); ++q) {
          index_type j = q.index(); // j >= i thanks to compact()
          value_type b_jk = *q;
          c(i,j) += b_ik*a_kk*b_jk;
        }
      }

      for (int l=k+1; l<n; ++l) {
        value_type a_kl = *a_++;
        for (const_row_iterator p = col(k).begin(); p != col(k).end(); ++p) {
          index_type i = p.index();
          value_type b_ik = *p;
          for (const_row_iterator q = col(l).begin(); q != col(l).end(); ++q) {
            index_type j = q.index();
            value_type b_jl = *q;
            value_type s = b_ik*a_kl*b_jl;
            if      (i < j) c(i,j) += s;
            else if (i > j) c(j,i) += s;
            else            c(i,j) += 2*s;
          }
        }
      }
    }
    return result;
  }

  /// A^T A where A is this matrix
  matrix this_transpose_times_this() const
  {
    index_type n = n_cols();
    matrix result(n, n);
    for (index_type i=0; i<n; ++i) {
      result(i, i) = col(i)*col(i);
      for (index_type j=i+1; j<n; ++j) {
        result(j, i) = result(i, j) = col(i)*col(j);
      }
    }
    result.compact();
    return result;
  }

  /// A^T W A where W is a diagonal matrix
  matrix
  this_transpose_times_diagonal_times_this(af::const_ref<value_type> const &w)
  const
  {
    SCITBX_ASSERT(w.size() == n_rows())(w.size())(n_rows());
    index_type n = n_cols();
    matrix result(n, n);
    for (index_type i=0; i<n; ++i) {
      result(i, i) = weighted_dot(col(i), w, col(i));
      for (index_type j=i+1; j<n; ++j) {
        result(j,i) = result(i,j) = weighted_dot(col(i), w, col(j));
      }
    }
    result.compact();
    return result;
  }

private:
  typedef af::shared<column_type> container_type;
  typedef af::ref<column_type> ref_type;

  typedef typename std::vector<index_type>::const_iterator const_row_idx_iter;
  index_type n_rows_;
  container_type columns;
  ref_type column;
};


template <class T>
inline
void matrix_times_dense_vector<T>::assign_to(af::ref<T> const &w) const {
  std::fill(w.begin(), w.end(), 0);
  for (int j=0; j < a.n_cols(); ++j) {
    for (typename matrix<T>::const_row_iterator p=a.col(j).begin();
         p != a.col(j).end();
         ++p)
    {
      typename matrix<T>::index_type i = p.index();
      typename matrix<T>::value_type a_ij = *p;
      w[i] += a_ij * v[j];
    }
  }
}

template <class T>
inline
void matrix_transpose_times_dense_vector<T>::assign_to(af::ref<T> const &w) const
{
  for (int j=0; j < a.n_cols(); ++j) {
    w[j] = 0;
    for (typename matrix<T>::const_row_iterator p=a.col(j).begin();
         p != a.col(j).end();
         ++p)
    {
      typename matrix<T>::index_type i = p.index();
      typename matrix<T>::value_type a_ij = *p;
      w[j] += a_ij * v[i];
    }
  }
}

}} // namespace scitbx::sparse


#endif
