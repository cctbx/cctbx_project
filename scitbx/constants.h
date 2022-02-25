/*! \file
    General purpose mathematical or physical %constants.
 */

#ifndef SCITBX_CONSTANTS_H
#define SCITBX_CONSTANTS_H

#include <cmath>

namespace scitbx {

  //! General purpose mathematical or physical %constants.
  namespace constants {
    //! mathematical constant pi
    static const double pi = 4. * std::atan(1.);
    //! mathematical constant pi*pi
    static const double pi_sq = pi * pi;
    //! mathematical constant 2*pi
    static const double two_pi = 8. * std::atan(1.);
    //! mathematical constant 2*pi*pi
    static const double two_pi_sq = 2. * pi_sq;
    //! mathematical constant 4*pi
    static const double four_pi = 16. * std::atan(1.);
    //! mathematical constant 4*pi*pi
    static const double four_pi_sq = two_pi * two_pi;
    //! mathematical constant 8*pi*pi
    static const double eight_pi_sq = 2. * four_pi_sq;
    //! mathematical constant pi/2
    static const double pi_2 = 2. * std::atan(1.);
    //! mathematical constant pi/180
    static const double pi_180 = std::atan(1.) / 45.;
    //! mathematical constant the Planck constant
    static const double Planck = 6.62607015e-34;
    //! mathematical constant neutron mass
    static const double neutron_mass = 1.67492749804e-27;
    //! mathematical constant neutron mass
    static const double m_n = 1.67492749804e-27;

    //! Factor for keV <-> Angstrom conversion.
    /*!
      https://www.bipm.org/utils/common/pdf/CGPM-2018/26th-CGPM-Resolutions.pdf

      h = Plank's Constant = 6.62607015e-34 J s
      c = speed of light = 2.99792458e+8 m/s
      1 keV = 1.e+3 * 1.602176634e-19 J
      1 A = Angstrom = 1.e-10 m

      E = (h * c) / lamda;

      Exponents: (-34 + 8) - (3 - 19 - 10) = 0
     */
    static const double
    factor_kev_angstrom = 6.62607015 * 2.99792458 / 1.602176634;
    static const double
    factor_ev_angstrom  = factor_kev_angstrom * 1000;
  }

  //! Conversions from degrees to radians.
  inline double deg_as_rad(double deg) { return deg * constants::pi_180;}
  //! Conversions from radians to degrees.
  inline double rad_as_deg(double rad) { return rad / constants::pi_180;}

} // namespace scitbx

#endif // SCITBX_CONSTANTS_H
