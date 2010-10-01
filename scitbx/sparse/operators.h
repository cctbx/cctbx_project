#ifndef SCITBX_SPARSE_OPERATORS_H
#define SCITBX_SPARSE_OPERATORS_H

namespace scitbx { namespace sparse { namespace operators {

  /// To implement af::ref<T, AccessorType> = sparse::matrix or sparse::vector
  /** Assignment or construction */
  template <class T>
  struct assign
  {
    void init(T *first, T *last) const {
      std::fill(first, last, T(0));
    }

    void operator()(T &b_ij, T a_ij) const { b_ij = a_ij; }
  };

  /// To implement af::ref<T, AccessorType> += sparse::matrix or sparse::vector
  template <class T>
  struct plus_equal
  {
    void init(T *first, T *last) const {};
    void operator()(T &b_ij, T a_ij) const { b_ij += a_ij; }
  };

  /// To implement af::ref<T, AccessorType> -= sparse::matrix or sparse::vector
  template <class T>
  struct minus_equal
  {
    void init(T *first, T *last) const {};
    void operator()(T &b_ij, T a_ij) const { b_ij -= a_ij; }
  };

  /// Whole matrix selection
  struct whole
  {
    bool operator()(int i, int j) const { return true; }
  };

  /// Upper diagonal selection
  struct upper_diagonal
  {
    bool operator()(int i, int j) const { return i <= j; }
  };

  /// Sparse matrix operating on a dense matrix
  template <class T, class IndexSelection, class Heir>
  class matrix_operating_on_dense_matrix
  {
  protected:
    template <class Operator, class AccessorType>
    void operate_on(Operator const &op, af::ref<T, AccessorType> const &b) const
    {
      Heir const &a = static_cast<Heir const &> (*this);
      SCITBX_ASSERT(a.n_cols() == b.n_columns() && a.n_rows() == b.n_rows())
                   (a.n_cols())(b.n_columns())(a.n_rows())(b.n_rows());
      IndexSelection select;
      op.init(b.begin(), b.end());
      a.compact();
      for (int j=0; j<a.n_cols(); ++j) {
        for (typename Heir::const_row_iterator p=a.col(j).begin();
             p != a.col(j).end();
             ++p)
        {
          typename Heir::index_type i = p.index();
          if (select(i, j)) {
            typename Heir::value_type a_ij = *p;
            op(b(i,j), a_ij);
          }
        }
      }
    }
  };

}}}

#endif // GUARD
