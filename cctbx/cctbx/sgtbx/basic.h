// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_BASIC_H
#define CCTBX_SGTBX_BASIC_H

#include <exception>
#include <string>
#include <cctbx/error.h>
#include <cctbx/basic/matrixlite.h>

//! Space Group Toolbox namespace.
namespace sgtbx {

  using namespace cctbx;

  class error_base_factor : public error {
    public:
      error_base_factor(const char* file, long line,
                        const std::string& msg = "") throw()
        : error(file, line, msg) {}
      virtual ~error_base_factor() throw() {}
  };

  using MatrixLite::itype::Mx33;
  using MatrixLite::itype::Vec3;

  // Maximum number of representative rotation matrices for
  // 3-dimensional crystallographic space groups.
  static const int nMaxReprRotMx = 24;

  const int STBF = 12; // Seitz           Matrix Translation Base Factor
  const int CRBF = 12; // Change of Basis Matrix Rotation    Base Factor
  const int CTBF = 72; // Change of Basis Matrix Translation Base Factor

  inline void sanity_check() {
    cctbx_assert(STBF % 12 == 0);
    cctbx_assert(CTBF >= 2 * STBF);
    cctbx_assert(CTBF % STBF == 0);
  }

} // namespace sgtbx

#endif // CCTBX_SGTBX_BASIC_H
