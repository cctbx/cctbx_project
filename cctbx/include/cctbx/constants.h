// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    General purpose mathematical or physical %constants.
 */

#ifndef CCTBX_CONSTANTS_H
#define CCTBX_CONSTANTS_H

#include <cctbx/fixes/cmath>

namespace cctbx {

  //! General purpose mathematical or physical %constants.
  namespace constants {
    //! mathematical constant pi
    static const double pi     = 4. * std::atan(1.);
    //! mathematical constant pi/2
    static const double pi_2   = 2. * std::atan(1.);
    //! mathematical constant pi/180
    static const double pi_180 = std::atan(1.) / 45.;

    //! Factor for keV <-> Angstrom conversion.
    /*!
      http://physics.nist.gov/PhysRefData/codata86/table2.html

      h = Plank's Constant = 6.6260755e-34 J s
      c = speed of light = 2.99792458e+8 m/s
      1 keV = 1.e+3 * 1.60217733e-19 J
      1 A = Angstrom = 1.e-10 m

      E = (h * c) / lamda;

      Exponents: (-34 + 8) - (3 - 19 - 10) = 0
     */
    static const double
    factor_keV_Angstrom = 6.6260755 * 2.99792458 / 1.60217733;
  }

} // namespace cctbx

#endif // CCTBX_CONSTANTS_H
