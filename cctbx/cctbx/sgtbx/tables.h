// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_TABLES_H
#define CCTBX_SGTBX_TABLES_H

#include <cctbx/sgtbx/matrix.h>

namespace sgtbx {
  namespace tables {

    namespace RotationMatrices
    {
      static const RotMx R_1_000(
         1,  0,  0,
         0,  1,  0,
         0,  0,  1 );
      static const RotMx R_2_001(
        -1,  0,  0,
         0, -1,  0,
         0,  0,  1 );
      static const RotMx R_2_1b0(
         0, -1,  0,
        -1,  0,  0,
         0,  0, -1 );
      static const RotMx R_2_110(
         0,  1,  0,
         1,  0,  0,
         0,  0, -1 );
      static const RotMx R_3_001(
         0, -1,  0,
         1, -1,  0,
         0,  0,  1 );
      static const RotMx R_3_111(
         0,  0,  1,
         1,  0,  0,
         0,  1,  0 );
      static const RotMx R_3i111(
         0,  1,  0,
         0,  0,  1,
         1,  0,  0 );
      static const RotMx R_4_001(
         0, -1,  0,
         1,  0,  0,
         0,  0,  1 );
      static const RotMx R_4i001(
         0,  1,  0,
        -1,  0,  0,
         0,  0,  1 );
      static const RotMx R_6_001(
         1, -1,  0,
         1,  0,  0,
         0,  0,  1 );
    }

  } // namespace tables
} // namespace sgtbx

#endif // CCTBX_SGTBX_TABLES_H
