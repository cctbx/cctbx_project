#ifndef SCITBX_MATRIX_MOVE_H
#define SCITBX_MATRIX_MOVE_H

#include <scitbx/matrix/packed.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace scitbx { namespace matrix {

  template <typename T>
  af::shared<T>
  copy_column(
    af::const_ref<T, af::c_grid<2> > const& self,
    unsigned i_column)
  {
    unsigned n_rows = self.accessor()[0];
    unsigned n_columns = self.accessor()[1];
    SCITBX_ASSERT(i_column < n_columns);
    af::shared<T> result(n_rows, af::init_functor_null<T>());
    T* r = result.begin();
    const T* s = &self[i_column];
    for(;n_rows>0;n_rows--) {
      *r++ = *s;
      s += n_columns;
    }
    return result;
  }

  template <typename T>
  af::versa<T, af::c_grid<2> >
  copy_block(
    af::const_ref<T, af::c_grid<2> > const& self,
    unsigned i_row,
    unsigned i_column,
    unsigned n_rows,
    unsigned n_columns)
  {
    unsigned self_n_rows = self.accessor()[0];
    unsigned self_n_columns = self.accessor()[1];
    SCITBX_ASSERT(i_row + n_rows <= self_n_rows);
    SCITBX_ASSERT(i_column + n_columns <= self_n_columns);
    af::versa<T, af::c_grid<2> >
      result(
        af::c_grid<2>(n_rows, n_columns),
        af::init_functor_null<T>());
    T* r = result.begin();
    const T* s = &self[i_row * self_n_columns + i_column];
    for(;n_rows>0;n_rows--) {
      r = std::copy(s, s+n_columns, r);
      s += self_n_columns;
    }
    return result;
  }

  template <typename T>
  void
  paste_block_in_place(
    af::ref<T, af::c_grid<2> > const& self,
    af::const_ref<T, af::c_grid<2> > const& block,
    unsigned i_row,
    unsigned i_column)
  {
    unsigned self_n_rows = self.accessor()[0];
    unsigned self_n_columns = self.accessor()[1];
    unsigned block_n_rows = block.accessor()[0];
    unsigned block_n_columns = block.accessor()[1];
    SCITBX_ASSERT(i_row + block_n_rows <= self_n_rows);
    SCITBX_ASSERT(i_column + block_n_columns <= self_n_columns);
    const T* b = block.begin();
    T* s = &self[i_row * self_n_columns + i_column];
    for(;block_n_rows>0;block_n_rows--) {
      std::copy(b, b+block_n_columns, s);
      b += block_n_columns;
      s += self_n_columns;
    }
  }

  template <typename T>
  void paste_column_in_place(af::ref<T, af::c_grid<2> > const &self,
                             af::const_ref<T> const &col,
                             unsigned j)
  {
    SCITBX_ASSERT(self.n_rows() == col.size())(self.n_rows())(col.size());
    SCITBX_ASSERT(j < self.n_columns())(j);
    for (unsigned i=0; i<self.n_rows(); ++i) self(i, j) = col(i);
  }

  template <typename T>
  void
  copy_upper_to_lower_triangle_in_place(
    af::ref<T, af::c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    unsigned n = a.accessor()[0];
    unsigned nn = n*n;
    unsigned ij = 0;
    unsigned ji0 = n;
    for(unsigned i=1;i<n;i++,ji0+=n) {
      ij += i;
      for(unsigned ji=ji0++;ji<nn;ji+=n) {
        a[ji] = a[ij++];
      }
    }
  }

  template <typename T>
  void
  copy_lower_to_upper_triangle_in_place(
    af::ref<T, af::c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    unsigned n = a.accessor()[0];
    unsigned i0 = 0;
    for(unsigned i=0;i<n;i++,i0+=n) {
      unsigned ij = i0;
      unsigned ji = i;
      for(unsigned j=0;j<i;j++,ji+=n) {
        a[ji] = a[ij++];
      }
    }
  }

  /// A zero matrix whose upper diagonal part is that of a.
  /** a must have more rows than columns */
  template <typename T>
  inline
  af::versa<T, af::mat_grid>
  copy_upper_triangle(af::const_ref<T, af::mat_grid> const &a) {
    int m=a.n_rows(), n=a.n_columns();
    SCITBX_ASSERT(m >= n);
    af::versa<T, af::mat_grid> result(af::mat_grid(n, n),
                                      af::init_functor_null<T>());
    for (int i=0; i<n; ++i) {
      std::fill_n(&result(i,0), i, T(0));
      std::copy(&a(i,i), &a(i,n), &result(i,i));
    }
    return result;
  }

  /// A zero matrix whose lower diagonal part is that of a.
  /** a must have less rows than columns */
  template <typename T>
  inline
  af::versa<T, af::mat_grid>
  copy_lower_triangle(af::const_ref<T, af::mat_grid> const &a) {
    int m=a.n_rows(), n=a.n_columns();
    SCITBX_ASSERT(m <= n);
    af::versa<T, af::mat_grid> result(af::mat_grid(m, m),
                                      af::init_functor_null<T>());
    for (int i=0; i<m; ++i) {
      std::fill(&result(i,i+1), &result(i,m), T(0));
      std::copy(&a(i,0), &a(i,i+1), &result(i,0));
    }
    return result;
  }

  template <typename T>
  void
  swap_rows_in_place(
    af::ref<T, af::c_grid<2> > const& a,
    unsigned i,
    unsigned j)
  {
    unsigned nr = a.accessor()[0];
    unsigned nc = a.accessor()[1];
    SCITBX_ASSERT(i < nr);
    SCITBX_ASSERT(j < nr);
    if (i == j) return;
    unsigned ik = i*nc;
    unsigned jk = j*nc;
    for(unsigned k=0;k<nc;k++) {
      std::swap(a[ik++], a[jk++]);
    }
  }

  template <typename T>
  void
  swap_columns_in_place(
    af::ref<T, af::c_grid<2> > const& a,
    unsigned i,
    unsigned j)
  {
    unsigned nc = a.accessor()[1];
    unsigned nrnc = a.accessor()[0]*nc;
    SCITBX_ASSERT(i < nc);
    SCITBX_ASSERT(j < nc);
    if (i == j) return;
    unsigned kj = j;
    for(unsigned ki=i;ki<nrnc;ki+=nc,kj+=nc) {
      std::swap(a[ki], a[kj]);
    }
  }

  template <typename T>
  void
  symmetric_upper_triangle_swap_rows_and_columns_in_place(
    af::ref<T, af::c_grid<2> > const& a,
    unsigned i,
    unsigned j)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    unsigned n = a.accessor()[0];
    SCITBX_ASSERT(i < n);
    SCITBX_ASSERT(j < n);
    if (i == j) return;
    if (i > j) std::swap(i, j);
    unsigned ki = i;
    unsigned kj = j;
    unsigned k = 0;
    for(;k<i;k++,ki+=n,kj+=n) {
      std::swap(a[ki], a[kj]);
    }
    unsigned ik = ki+1;
    for(k++;k<j;k++) {
      kj += n;
      std::swap(a[ik++], a[kj]);
    }
    kj += n;
    std::swap(a[ki], a[kj]);
    for(k++;k<n;k++) {
      std::swap(a[++ik], a[++kj]);
    }
  }

  template <typename T>
  void
  packed_u_swap_rows_and_columns_in_place(
    af::ref<T> const& u,
    unsigned i,
    unsigned j)
  {
    unsigned n = af::dimension_from_packed_size(u.size());
    SCITBX_ASSERT(i < n);
    SCITBX_ASSERT(j < n);
    if (i == j) return;
    if (i > j) std::swap(i, j);
    unsigned d = j - i;
    unsigned ki = i;
    unsigned k=0;
    while (k<i) {
      std::swap(u[ki], u[ki+d]);
      ki += n-(++k);
    }
    unsigned ik = ki;
    unsigned kj = ki + n-(++k) + d;
    while (k<j) {
      std::swap(u[++ik], u[kj]);
      kj += n-(++k);
    }
    std::swap(u[ki], u[kj]);
    ik++;
    for(k++;k<n;k++) {
      std::swap(u[++ik], u[++kj]);
    }
  }

  /// Copy the diagonal and superdiagonal of matrix src into matrix dst
  template <typename T1, class A1, typename T2, class A2>
  void copy_upper_bidiagonal(af::ref<T1, A1>       const &dst,
                             af::const_ref<T2, A2> const &src)
  {
    int p = std::min(src.n_rows(), src.n_columns());
    for (int i=0; i<p; ++i){
      dst(i,i) = src(i,i);
      if (i < p-1) dst(i,i+1) = src(i,i+1);
    }
  }

  /// Copy the diagonal and subdiagonal of matrix src into matrix dst
  template <typename T1, class A1, typename T2, class A2>
  void copy_lower_bidiagonal(af::ref<T1, A1>       const &dst,
                             af::const_ref<T2, A2> const &src)
  {
    int p = std::min(src.n_rows(), src.n_columns());
    for (int i=0; i<p; ++i){
      dst(i,i) = src(i,i);
      if (i < p-1) dst(i+1,i) = src(i+1,i);
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_MOVE_H
