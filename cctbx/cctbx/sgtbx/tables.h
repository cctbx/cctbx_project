// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
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

    namespace CrystalSystem {

      enum Code {
        Undefined,
        Triclinic,
        Monoclinic,
        Orthorhombic,
        Tetragonal,
        Trigonal,
        Hexagonal,
        Cubic
      };

      inline const char* LabelOf(const Code& c) {
        static const char* Names[] = {
          "Undefined",
          "Triclinic",
          "Monoclinic",
          "Orthorhombic",
          "Tetragonal",
          "Trigonal",
          "Hexagonal",
          "Cubic",
        };
        return Names[c];
      }

    } // namespace CrystalSystem

    namespace MatrixGroup {

      class Code {
        public:
          inline Code(const CrystalSystem::Code& X, int L, int P, int M)
            : m_X(X), m_L(L), m_P(P), m_M(M) {}
          inline bool operator==(const Code& rhs) const {
            return    m_X == rhs.m_X
                   && m_L == rhs.m_L
                   && m_P == rhs.m_P
                   && m_M == rhs.m_M;
          }
          inline bool operator!=(const Code& rhs) const {
            return !(*this == rhs);
          }
          inline int index() const { return m_M; }
          inline Code PointGroupType() const {
            return Code(m_X, m_L - m_P, 0, m_M + m_P);
          }
          inline Code LaueGroupType() const {
            return Code(m_X, 0, 0, m_M + m_L);
          }
          inline CrystalSystem::Code CrystalSystem() const {
            return m_X;
          }
          inline const char* Label() const {
            static const char* Names[] = {
              "Undefined",
              "Unknown",
              "1",
              "-1",
              "2",
              "m",
              "2/m",
              "222",
              "mm2",
              "mmm",
              "4",
              "-4",
              "4/m",
              "422",
              "4mm",
              "-4m2",
              "-42m",
              "4/mmm",
              "3",
              "-3",
              "321",
              "312",
              "32",
              "3m1",
              "31m",
              "3m",
              "-3m1",
              "-31m",
              "-3m",
              "6",
              "-6",
              "6/m",
              "622",
              "6mm",
              "-6m2",
              "-62m",
              "6/mmm",
              "23",
              "m-3",
              "432",
              "-43m",
              "m-3m"
            };
            return Names[m_M];
          }
        private:
          CrystalSystem::Code m_X;
          int m_L;
          int m_P;
          int m_M;
      };

      static const Code Undefined  (CrystalSystem::Undefined,     0,  0,  0);
      static const Code Unknown    (CrystalSystem::Undefined,     0,  0,  1);
      static const Code MGC_1      (CrystalSystem::Triclinic,     1,  0,  2);
      static const Code MGC_1b     (CrystalSystem::Triclinic,     0,  0,  3);
      static const Code MGC_2      (CrystalSystem::Monoclinic,    2,  0,  4);
      static const Code MGC_m      (CrystalSystem::Monoclinic,    1,  0,  5);
      static const Code MGC_2_m    (CrystalSystem::Monoclinic,    0,  0,  6);
      static const Code MGC_222    (CrystalSystem::Orthorhombic,  2,  0,  7);
      static const Code MGC_mm2    (CrystalSystem::Orthorhombic,  1,  0,  8);
      static const Code MGC_mmm    (CrystalSystem::Orthorhombic,  0,  0,  9);
      static const Code MGC_4      (CrystalSystem::Tetragonal,    2,  0, 10);
      static const Code MGC_4b     (CrystalSystem::Tetragonal,    1,  0, 11);
      static const Code MGC_4_m    (CrystalSystem::Tetragonal,    0,  0, 12);
      static const Code MGC_422    (CrystalSystem::Tetragonal,    4,  0, 13);
      static const Code MGC_4mm    (CrystalSystem::Tetragonal,    3,  0, 14);
      static const Code MGC_4b2m   (CrystalSystem::Tetragonal,    2,  1, 15);
      static const Code MGC_4bm2   (CrystalSystem::Tetragonal,    1,  0, 16);
      static const Code MGC_4_mmm  (CrystalSystem::Tetragonal,    0,  0, 17);
      static const Code MGC_3      (CrystalSystem::Trigonal,      1,  0, 18);
      static const Code MGC_3b     (CrystalSystem::Trigonal,      0,  0, 19);
      static const Code MGC_321    (CrystalSystem::Trigonal,      8,  2, 20);
      static const Code MGC_312    (CrystalSystem::Trigonal,      7,  1, 21);
      static const Code MGC_32     (CrystalSystem::Trigonal,      6,  0, 22);
      static const Code MGC_3m1    (CrystalSystem::Trigonal,      5,  2, 23);
      static const Code MGC_31m    (CrystalSystem::Trigonal,      4,  1, 24);
      static const Code MGC_3m     (CrystalSystem::Trigonal,      3,  0, 25);
      static const Code MGC_3bm1   (CrystalSystem::Trigonal,      2,  2, 26);
      static const Code MGC_3b1m   (CrystalSystem::Trigonal,      1,  1, 27);
      static const Code MGC_3bm    (CrystalSystem::Trigonal,      0,  0, 28);
      static const Code MGC_6      (CrystalSystem::Hexagonal,     2,  0, 29);
      static const Code MGC_6b     (CrystalSystem::Hexagonal,     1,  0, 30);
      static const Code MGC_6_m    (CrystalSystem::Hexagonal,     0,  0, 31);
      static const Code MGC_622    (CrystalSystem::Hexagonal,     4,  0, 32);
      static const Code MGC_6mm    (CrystalSystem::Hexagonal,     3,  0, 33);
      static const Code MGC_6b2m   (CrystalSystem::Hexagonal,     2,  1, 34);
      static const Code MGC_6bm2   (CrystalSystem::Hexagonal,     1,  0, 35);
      static const Code MGC_6_mmm  (CrystalSystem::Hexagonal,     0,  0, 36);
      static const Code MGC_23     (CrystalSystem::Cubic,         1,  0, 37);
      static const Code MGC_m3b    (CrystalSystem::Cubic,         0,  0, 38);
      static const Code MGC_432    (CrystalSystem::Cubic,         2,  0, 39);
      static const Code MGC_4b3m   (CrystalSystem::Cubic,         1,  0, 40);
      static const Code MGC_m3bm   (CrystalSystem::Cubic,         0,  0, 41);

    } // namespace MatrixGroup
  } // namespace tables
} // namespace sgtbx

#endif // CCTBX_SGTBX_TABLES_H
