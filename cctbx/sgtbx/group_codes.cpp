/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/group_codes.h>

namespace cctbx { namespace sgtbx {

  namespace crystal_system {

    const char* label(code const& c)
    {
      static const char* names[] = {
        "Undefined",
        "Triclinic",
        "Monoclinic",
        "Orthorhombic",
        "Tetragonal",
        "Trigonal",
        "Hexagonal",
        "Cubic",
      };
      return names[c];
    }

  } // namespace crystal_system

  namespace matrix_group {

    const char* code::label() const
    {
      static const char* names[] = {
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
      return names[m_];
    }

  } // namespace matrix_group

}} // namespace cctbx::sgtbx
