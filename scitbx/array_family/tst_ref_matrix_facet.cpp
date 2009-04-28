#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/ref_matrix.h>
#include <iostream>

using namespace scitbx;

typedef af::c_grid<2> dim;
typedef af::const_ref<int, dim> mat_const_ref_t;
typedef af::ref<int, dim> mat_ref_t;

typedef af::tiny<int, 6>  a6_t;
typedef af::tiny<int, 8>  a8_t;
typedef af::tiny<int, 9>  a9_t;
typedef af::tiny<int, 12> a12_t;

int main() {
  {
    a6_t u(1,2,3,4,5,6);
    mat_const_ref_t a(u.begin(), 2,3);
    SCITBX_ASSERT(a.n_rows() == 2)(a.n_rows());
    SCITBX_ASSERT(a.n_columns() == 3)(a.n_columns());
    SCITBX_ASSERT(!a.is_square());
  }
  {
    a6_t u(1,2,3,4,5,6);
    mat_ref_t b(u.begin(), 2,3);
    b += 2;
    SCITBX_ASSERT((u.all_eq(a6_t(3,4,5,6,7,8))));
    b /=2;
    SCITBX_ASSERT((u.all_eq(a6_t(1,2,2,3,3,4))));
  }
  {
    a6_t ta(1,2,3,4,5,6);
    a12_t tb;
    for(int i=0;i<12;i++) tb[i] = i+1;
    a8_t tab;
    mat_ref_t a(ta.begin(), 2, 3);
    mat_ref_t b(tb.begin(), 3, 4);
    mat_ref_t ab(tab.begin(), 2, 4);
    af::multiply(a, b, ab);
    SCITBX_ASSERT(ab.all_eq(
        mat_const_ref_t(a8_t(38,44,50,56,83,98,113,128).begin(), 2, 4)));
  }
  {
    a6_t ta(1,2,3,4,5,6);
    a6_t tb(1,2,3,4);
    a6_t tc(2,0,0,2);
    mat_ref_t a(ta.begin(), 2, 3);
    mat_ref_t b(tb.begin(), 2, 2);
    mat_ref_t c(tc.begin(), 2, 2);
    SCITBX_ASSERT(!a.is_square());
    SCITBX_ASSERT(!a.is_diagonal());
    SCITBX_ASSERT(b.is_square());
    SCITBX_ASSERT(!b.is_diagonal());
    SCITBX_ASSERT(c.is_square());
    SCITBX_ASSERT(c.is_diagonal());
  }
  {
    a6_t ta(1,2,3,4,5,6);
    mat_ref_t a(ta.begin(), 3, 2);
    a.swap_rows(0, 1);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(3,4,1,2,5,6).begin(), 3, 2)));
    a.swap_rows(1, 2);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(3,4,5,6,1,2).begin(), 3, 2)));
    a.swap_rows(0, 2);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(1,2,5,6,3,4).begin(), 3, 2)));
    a.swap_rows(2, 1);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(1,2,3,4,5,6).begin(), 3, 2)));
    a.swap_columns(0, 1);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(2,1,4,3,6,5).begin(), 3, 2)));
    a.swap_columns(1, 0);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(1,2,3,4,5,6).begin(), 3, 2)));
    a = mat_ref_t(a.begin(), 2, 3);
    a.swap_columns(2, 1);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(1,3,2,4,6,5).begin(), 2, 3)));
    a.swap_columns(0, 2);
    SCITBX_ASSERT(a.all_eq(mat_const_ref_t(a6_t(2,3,1,5,6,4).begin(), 2, 3)));
  }
  {
    a9_t ta;
    mat_ref_t a(ta.begin(), 3, 3);
    a.set_diagonal(3);
    SCITBX_ASSERT(ta.all_eq(a9_t(3,0,0,0,3,0,0,0,3)));
    SCITBX_ASSERT(a.is_diagonal());
    a.set_identity();
    SCITBX_ASSERT(ta.all_eq(a9_t(1,0,0,0,1,0,0,0,1)));
    SCITBX_ASSERT(a.is_diagonal());
  }
  {
    a9_t ta(1,2,3,4,5,6,7,8,9);
    mat_ref_t a(ta.begin(), 3, 3);
    a.transpose_in_place();
    SCITBX_ASSERT(ta.all_eq(a9_t(1,4,7,2,5,8,3,6,9)));
    a.transpose_in_place();
    SCITBX_ASSERT(ta.all_eq(a9_t(1,2,3,4,5,6,7,8,9)));
  }
  {
    a6_t ta(1,2,3,4,5,6);
    mat_ref_t a(ta.begin(), 2, 3);
    a.transpose_in_place();
    SCITBX_ASSERT(ta.all_eq(a6_t(1,4,2,5,3,6)));
    a.transpose_in_place();
    SCITBX_ASSERT(ta.all_eq(a6_t(1,2,3,4,5,6)));
  }

  std::cout << "OK\n";
}
