#include <cctbx/sgtbx/tensor_rank_2.h>

namespace cctbx { namespace sgtbx { namespace tensor_rank_2 {

  int*
  constraints_raw(const int* r, int* c)
  {
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
    return c;
  }

  int*
  constraints_raw(
    af::const_ref<rt_mx> const& symmetry_matrices,
    std::size_t i_first_matrix_to_use,
    bool reciprocal_space,
    int* c)
  {
    rot_mx r_transpose;
    const int* r;
    if (reciprocal_space) r = r_transpose.num().begin();
    for(std::size_t i=i_first_matrix_to_use;i<symmetry_matrices.size();i++) {
      if (reciprocal_space) {
        r_transpose = symmetry_matrices[i].r().transpose();
      }
      else {
        r = symmetry_matrices[i].r().num().begin();
      }
      c = constraints_raw(r, c);
    }
    return c;
  }

}}} // namespace cctbx::sgtbx::tensor_rank_2
