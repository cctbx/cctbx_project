#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/mat_ref.h>

namespace cctbx { namespace sgtbx {

  af::versa<int, af::c_grid<2> >
  space_group
  ::tensor_constraints(bool reciprocal_space) const
  {
    af::c_grid<2> result_accessor((n_smx()-1)*6, 6);
    af::versa<int, af::c_grid<2> >
      result(result_accessor, af::init_functor_null<int>());
    rot_mx r_transpose;
    const int* r;
    if (reciprocal_space) r = r_transpose.num().begin();
    int* c = result.begin();
    for(std::size_t i_smx=1;i_smx<n_smx();i_smx++) {
      if (reciprocal_space) {
        r_transpose = smx_[i_smx].r().transpose();
      }
      else {
        r = smx_[i_smx].r().num().begin();
      }
      // See sgtbx/space_group.h for the Mathematica code used to
      // generate the following expressions.
      // row 0
      *c++ = r[0]*r[0]-1;
      *c++ = r[3]*r[3];
      *c++ = r[6]*r[6];
      *c++ = 2*r[0]*r[3];
      *c++ = 2*r[0]*r[6];
      *c++ = 2*r[3]*r[6];
      // row 1
      *c++ = r[1]*r[1];
      *c++ = r[4]*r[4]-1;
      *c++ = r[7]*r[7];
      *c++ = 2*r[1]*r[4];
      *c++ = 2*r[1]*r[7];
      *c++ = 2*r[4]*r[7];
      // row 2
      *c++ = r[2]*r[2];
      *c++ = r[5]*r[5];
      *c++ = r[8]*r[8]-1;
      *c++ = 2*r[2]*r[5];
      *c++ = 2*r[2]*r[8];
      *c++ = 2*r[5]*r[8];
      // row 3
      *c++ = r[0]*r[1];
      *c++ = r[3]*r[4];
      *c++ = r[6]*r[7];
      *c++ = r[1]*r[3]+r[0]*r[4]-1;
      *c++ = r[1]*r[6]+r[0]*r[7];
      *c++ = r[4]*r[6]+r[3]*r[7];
      // row 4
      *c++ = r[0]*r[2];
      *c++ = r[3]*r[5];
      *c++ = r[6]*r[8];
      *c++ = r[2]*r[3]+r[0]*r[5];
      *c++ = r[2]*r[6]+r[0]*r[8]-1;
      *c++ = r[5]*r[6]+r[3]*r[8];
      // row 5
      *c++ = r[1]*r[2];
      *c++ = r[4]*r[5];
      *c++ = r[7]*r[8];
      *c++ = r[2]*r[4]+r[1]*r[5];
      *c++ = r[2]*r[7]+r[1]*r[8];
      *c++ = r[5]*r[7]+r[4]*r[8]-1;
    }
    CCTBX_ASSERT(c == result.end());
    scitbx::mat_ref<int> re_mx(result.begin(), (n_smx()-1)*6, 6);
    CCTBX_ASSERT(scitbx::matrix::row_echelon::form(re_mx) <= 6);
    result_accessor = af::c_grid<2>(re_mx.n_rows(), 6);
    result.resize(result_accessor);
    return result;
  }

}} // namespace cctbx::sgtbx
