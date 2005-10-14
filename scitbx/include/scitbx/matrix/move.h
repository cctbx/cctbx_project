#ifndef SCITBX_MATRIX_MOVE_H
#define SCITBX_MATRIX_MOVE_H

#include <scitbx/matrix/packed.h>

namespace scitbx { namespace matrix {

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
    unsigned nn = n*n;
    unsigned i0 = 0;
    for(unsigned i=0;i<n;i++,i0+=n) {
      unsigned ij = i0;
      unsigned ji = i;
      for(unsigned j=0;j<i;j++,ji+=n) {
        a[ji] = a[ij++];
      }
    }
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
    unsigned n = symmetric_n_from_packed_size(u.size());
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

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_MOVE_H
