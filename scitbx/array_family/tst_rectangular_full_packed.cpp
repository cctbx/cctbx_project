#include <scitbx/array_family/accessors/rectangular_full_packed.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/simple_io.h>
#include <fast_linalg/lapacke.h>
#include <stdio.h>

namespace af = scitbx::af;

typedef af::versa<double, af::packed_u_accessor> upper_packed_t;
typedef af::versa<double, af::packed_l_accessor> lower_packed_t;

typedef af::versa<double, af::rectangular_full_packed_accessor<af::upper> >
        upper_rfp_t;
typedef af::versa<double, af::rectangular_full_packed_accessor<af::lower> >
        lower_rfp_t;

typedef af::const_ref<double, af::f_grid<2,int> > rect_t;

int main() {
  if (!fast_linalg::is_initialised()) {
    printf("Test of accessor for rectangular full packed format skipped\n");
    return 0;
  }
  bool failed = false;
  for(int n=0; n<10; n++) {
    // simulate upper packed by columns
    lower_packed_t u_p(n);
    double *p = u_p.begin();
    for(int j=0; j<n; j++) for(int i=0; i<=j; i++) *p++ = 10*i + j;
    upper_rfp_t u_rfp(n);

    // Build RFP
    fast_linalg::tpttf(fast_linalg::LAPACK_COL_MAJOR, 'N', 'U',
      n, u_p.begin(), u_rfp.begin());

    // Check accessor works
    bool this_failed = false;
    for(int j=0; j<n; j++) for(int i=0; i<=j; i++) {
      if(u_rfp(i,j) != 10*i+j) {
        this_failed = true;
        break;
      }
    }
    if(this_failed) {
      failed = true;
      af::f_grid<2,int> rfp(n&1 ? n       : n+1,
                            n&1 ? (n+1)/2 : n/2);
      rect_t u_f(u_rfp.begin(), rfp);
      printf("n=%d\n---\n", n);
      for(int i=0; i<rfp.n_rows(); i++) {
        for(int j=0; j<rfp.n_columns(); j++) printf("%02d ", int(u_f(i,j)));
        printf("\n");
      }
      printf("\n");
      for(int i=0; i<n; i++) {
        for(int j=0; j<i; j++) printf("   ");
        for(int j=i; j<n; j++) printf("%02d ", int(u_rfp(i,j)));
        printf("\n");
      }
      printf("\n\n");
    }
  }
  if(failed) {
    printf("Failed\n");
    return 1;
  }
  else {
    printf("OK\n");
    return 0;
  }
}
