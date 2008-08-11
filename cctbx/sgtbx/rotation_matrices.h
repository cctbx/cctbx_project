#ifndef CCTBX_SGTBX_ROTATION_MATRICES_H
#define CCTBX_SGTBX_ROTATION_MATRICES_H

#include <cctbx/sgtbx/rot_mx.h>

namespace cctbx { namespace sgtbx { namespace tables {

  namespace rotation_matrices
  {
    static const rot_mx r_1_000(
       1,  0,  0,
       0,  1,  0,
       0,  0,  1 );
    static const rot_mx r_2_001(
      -1,  0,  0,
       0, -1,  0,
       0,  0,  1 );
    static const rot_mx r_2_1b0(
       0, -1,  0,
      -1,  0,  0,
       0,  0, -1 );
    static const rot_mx r_2_110(
       0,  1,  0,
       1,  0,  0,
       0,  0, -1 );
    static const rot_mx r_3_001(
       0, -1,  0,
       1, -1,  0,
       0,  0,  1 );
    static const rot_mx r_3_111(
       0,  0,  1,
       1,  0,  0,
       0,  1,  0 );
    static const rot_mx r_3i111(
       0,  1,  0,
       0,  0,  1,
       1,  0,  0 );
    static const rot_mx r_4_001(
       0, -1,  0,
       1,  0,  0,
       0,  0,  1 );
    static const rot_mx r_4i001(
       0,  1,  0,
      -1,  0,  0,
       0,  0,  1 );
    static const rot_mx r_6_001(
       1, -1,  0,
       1,  0,  0,
       0,  0,  1 );
  }

}}} // namespace cctbx::sgtbx::tables

#endif // CCTBX_SGTBX_ROTATION_MATRICES_H
