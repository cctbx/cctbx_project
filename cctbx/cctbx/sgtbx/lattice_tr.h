// $Id$

#ifndef CCTBX_SGTBX_LATTICE_TR_H
#define CCTBX_SGTBX_LATTICE_TR_H

#include <cctbx/sgtbx/matrix.h>

namespace sgtbx {
  namespace lattice {
    namespace tables {

      namespace ConventionalCentringTypes
      {
        static const TrVec P[] = { TrVec12(0, 0, 0),
                                 };
        static const TrVec A[] = { TrVec12(0, 0, 0),
                                   TrVec12(0, 6, 6),
                                };
        static const TrVec B[] = { TrVec12(0, 0, 0),
                                   TrVec12(6, 0, 6),
                                 };
        static const TrVec C[] = { TrVec12(0, 0, 0),
                                   TrVec12(6, 6, 0),
                                 };
        static const TrVec I[] = { TrVec12(0, 0, 0),
                                   TrVec12(6, 6, 6),
                                 };
        static const TrVec R[] = { TrVec12(0, 0, 0),
                                   TrVec12(8, 4, 4),
                                   TrVec12(4, 8, 8),
                                 };
        static const TrVec Q[] = { TrVec12(0, 0, 0), // reverse setting
                                   TrVec12(4, 8, 4), // for internal use
                                   TrVec12(8, 4, 8), // only
                                 };
        static const TrVec H[] = { TrVec12(0, 0, 0),
                                   TrVec12(8, 4, 0),
                                   TrVec12(4, 8, 0),
                                 };
        static const TrVec F[] = { TrVec12(0, 0, 0),
                                   TrVec12(0, 6, 6),
                                   TrVec12(6, 0, 6),
                                   TrVec12(6, 6, 0),
                                 };
      }

      namespace ConventionalZ2PMatrices // change-of-basis matrices
      {                                 // (coordinate transformations)
        static const RotMx PP(
           1,  0,  0,
           0,  1,  0,
           0,  0,  1 );
        static const RotMx AP(
          -1,  0,  0,
           0, -1,  1,
           0,  1,  1 );
        static const RotMx BP(
          -1,  0,  1,
           0, -1,  0,
           1,  0,  1 );
        static const RotMx CP(
           1,  1,  0,
           1, -1,  0,
           0,  0, -1 );
        static const RotMx IP(
           0,  1,  1,
           1,  0,  1,
           1,  1,  0 );
        static const RotMx RP(
           1,  0,  1,
          -1,  1,  1,
           0, -1,  1 );
        static const RotMx HP(
           1,  1,  0,
          -1,  2,  0,
           0,  0,  1 );
        static const RotMx FP(
          -1,  1,  1,
           1, -1,  1,
           1,  1, -1 );
      }

    } // namespace tables

    struct CentringTypeMap {
      char        Symbol;
      int         nTrs;
      const TrVec* Trs;
    };

    namespace tables {

      static const CentringTypeMap ConventionalCentringTypeMap[] = {
        { 'P', 1, ConventionalCentringTypes::P },
        { 'A', 2, ConventionalCentringTypes::A },
        { 'B', 2, ConventionalCentringTypes::B },
        { 'C', 2, ConventionalCentringTypes::C },
        { 'I', 2, ConventionalCentringTypes::I },
        { 'R', 3, ConventionalCentringTypes::R },
        { 'Q', 3, ConventionalCentringTypes::Q }, // reverse setting
        { 'H', 3, ConventionalCentringTypes::H },
        { 'F', 4, ConventionalCentringTypes::F },
        { '\0', 0, 0 },
      };

    } // namespace tables

    const CentringTypeMap* getConventionalCentringType(char Symbol);

  } // namespace lattice
} // namespace sgtbx

#endif // CCTBX_SGTBX_LATTICE_TR_H
