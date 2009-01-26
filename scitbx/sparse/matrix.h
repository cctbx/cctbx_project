#ifndef SCITBX_SPARSE_MATRIX_H
#define SCITBX_SPARSE_MATRIX_H

#include <algorithm>
#include <functional>
#include <vector>
#include <boost/lambda/bind.hpp>
#include <scitbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/sparse/vector.h>

namespace scitbx { namespace sparse {

template<class MatrixType, class VectorType>
struct matrix_x_vector
{
  typedef typename MatrixType::value_type value_type;
  typedef typename MatrixType::const_row_iterator const_row_iterator;
  typedef typename MatrixType::row_index row_index;
  typedef typename MatrixType::column_index column_index;

  std::vector<value_type> w;

  matrix_x_vector(row_index n_rows) : w(n_rows, 0) {}

  void operator()(MatrixType const &m,
                  VectorType const &v,
                  VectorType &u)
  {
    std::vector<row_index> nz;
    typedef typename std::vector<row_index>::const_iterator
            const_row_idx_iter;
    for (const_row_iterator pv=v.begin(); pv != v.end(); pv++) {
      column_index j = pv.index();
      value_type v_j = *pv;
      for (const_row_iterator pm=m.col(j).begin(); pm != m.col(j).end(); pm++) {
        row_index i = pm.index();
        value_type m_ij = *pm;
        if (w[i] == 0) nz.push_back(i);
        w[i] += m_ij * v_j;
      }
    }
    for (const_row_idx_iter p = nz.begin(); p != nz.end(); p++) {
      u[*p] = w[*p];
      w[*p] = 0; // ready to perform another product
    }
  }
};


/// A sparse matrix, represented by a sequence of sparse columns
/** All linear operations are therefore performed using column version of the
relevant algorithms, taking great care of never touching structurally zero
elements.

Also, those columns are kept unsorted (c.f. class sparse::vector) and the
precondition about duplicate applies here too.

When constructed with a definite number of rows, this number is retained and
immutable. However it is up to the user to make sure that no element is assigned
with an 1st index greater or equal to that number. Calling sort_indices will
however prune those illegal elements.

When constructed with an undefinite number of rows, that number stays so until
the first time the member function sort_indices is called, after which it is set
to the greatest size of all columns and it is immutable from then on.
This is to make it easy to fill a sparse matrix in a context
where std::vector::push_back or the like would normally be used to fill columns.
*/
template<class T>
class matrix
{
  private:
    /* Let's take advantage of scitbx reference counted array semantic...
    C.f., for example, member function "transpose"
    */
    typedef af::shared< vector<T> > container_type;

  public:
    typedef T value_type;
    typedef vector<T> column_type;
    typedef typename vector<T>::index_type row_index;
    typedef typename container_type::size_type column_index;
    typedef typename vector<T>::iterator row_iterator;
    typedef typename vector<T>::const_iterator const_row_iterator;

  public:
    /// Construct a zero matrix with the given number of rows and columns
    matrix(boost::optional<row_index> rows, column_index cols)
      : n_rows_(rows), column(cols)
    {
      for (column_index j=0; j < cols; j++) column[j] = column_type(rows);
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
    const typename vector<T>::element_reference
    operator()(row_index i, column_index j) const {
      return column[j][i];
    }

    /// Number of columns
    column_index n_cols() const {
      return column.size();
    }

    /// Number of rows
    row_index n_rows() const {
      SCITBX_ASSERT(n_rows_);
      return *n_rows_;
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
    /** C.f. the member function of same name in class scitbx::sparse::vector */
    void sort_indices() {
      using namespace boost::lambda;
      for (column_index j=0; j < n_cols(); j++) col(j).sort_indices();
      if (!n_rows_) {
        n_rows_ = std::max_element(
          column.begin(), column.end(),
          bind(&vector<T>::size, _1) < bind(&vector<T>::size, _2))->size();
      }
    }

    /// Permute the rows, in place
    template<class PermutationType>
    matrix& permute_rows(PermutationType const& permutation) {
      SCITBX_ASSERT(n_rows() == permutation.size())
                   ( n_rows() )( permutation.size() );
      for (column_index j=0; j < n_cols(); j++) col(j).permute(permutation);
      return *this;
    }

    /// Matrix times vector
    vector<T> operator*(vector<T> const& v) const {
      SCITBX_ASSERT(n_cols() == v.size())
                   ( n_cols() )( v.size() );
      vector<T> result(n_rows());
      matrix_x_vector<matrix, vector<T> > multiply(n_rows());
      multiply(*this, v, result);
      return result;
    }

    /// Dense row vector times matrix
    typedef typename vector<T>::dense_vector_const_ref dense_vector_const_ref;
    typedef af::shared<T> dense_vector;

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
      matrix result(a.n_rows(), b.n_cols());
      matrix_x_vector<matrix, vector<T> > multiply(a.n_rows());
      for (column_index j=0; j < b.n_cols(); j++) {
        multiply(a, b.col(j), result.col(j));
      }
      return result;
    }

  private:
    typedef typename std::vector<row_index>::const_iterator const_row_idx_iter;
    container_type column;
    mutable boost::optional<row_index> n_rows_;
};

}} // namespace scitbx::sparse


#endif
