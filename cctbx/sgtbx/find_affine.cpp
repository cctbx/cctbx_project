#include <cctbx/sgtbx/find_affine.h>

namespace cctbx { namespace sgtbx {

  find_affine::find_affine(space_group const& group, int range)
  {
    space_group tidy_group(group);
    tidy_group.make_tidy();
    // 0 1 2
    // 3 4 5
    // 6 7 8
    // Manually optimized search for matrices with determinant 1.
    for(int i4=-range;i4<=range;i4++) {
    for(int i8=-range;i8<=range;i8++) { int i48 = i4*i8;
    for(int i5=-range;i5<=range;i5++) {
    for(int i7=-range;i7<=range;i7++) { int i4857 = i48-i5*i7;
    for(int i3=-range;i3<=range;i3++) { int i38 = i3*i8;
                                        int i37 = i3*i7;
    for(int i6=-range;i6<=range;i6++) { int i3746 = i37-i4*i6;
                                        int i3856 = i38-i5*i6;
      if (i3746 == 0) {
        if (i3856 == 0) {
          if (i4857 == -1 || i4857 == 1) {
            for(int i1=-range;i1<=range;i1++) {
              for(int i2=-range;i2<=range;i2++) {
                rt_mx c(rot_mx(i4857,i1,i2,i3,i4,i5,i6,i7,i8));
                change_of_basis_op cb_op(c);
                if (   tidy_group.n_smx() == 1
                    || tidy_group.change_basis(cb_op) == tidy_group) {
                  cb_mx_.push_back(cb_op.c());
                }
              }
            }
          }
        }
        else {
          for(int i0=-range;i0<=range;i0++) {
            int i04857_one = i0*i4857 - 1;
            int i1 = i04857_one / i3856;
            if (-range <= i1 && i1 <= range && i1*i3856 == i04857_one) {
              for(int i2=-range;i2<=range;i2++) {
                rt_mx c(rot_mx(i0,i1,i2,i3,i4,i5,i6,i7,i8));
                change_of_basis_op cb_op(c);
                if (   tidy_group.n_smx() == 1
                    || tidy_group.change_basis(cb_op) == tidy_group) {
                  cb_mx_.push_back(cb_op.c());
                }
              }
            }
          }
        }
      }
      else {
        for(int i0=-range;i0<=range;i0++) {
          int i04857 = i0*i4857;
          for(int i1=-range;i1<=range;i1++) {
            int one_i04857_13856 = 1 - i04857 + i1*i3856;
            int i2 = one_i04857_13856 / i3746;
            if (-range <= i2 && i2 <= range && i2*i3746 == one_i04857_13856) {
              rt_mx c(rot_mx(i0,i1,i2,i3,i4,i5,i6,i7,i8));
              change_of_basis_op cb_op(c);
              if (   tidy_group.n_smx() == 1
                  || tidy_group.change_basis(cb_op) == tidy_group) {
                cb_mx_.push_back(cb_op.c());
              }
            }
          }
        }
      }
    }}}}}}
  }

}} // namespace cctbx::sgtbx
