// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#ifndef CCTBX_ELTBX_HENKE_H
#define CCTBX_ELTBX_HENKE_H

#include <string>
#include <cctbx/eltbx/efpfdp.h>

namespace eltbx {

  //! Access to Henke tables.
  /*! Henke tables are available for elements with Z=1-92.
      Each table contains 500+ points on a uniform logarithmic mesh
      from 10 to 30,000 eV with points added 0.1 eV above and below
      "sharp" absorption edges. The atomic scattering factors are
      based upon experimental measurements of the atomic
      photoabsorption cross section. The absorption measurements
      provide values for the imaginary part of the atomic scattering
      factor. The real part is calculated from the absorption
      measurements using the Kramers-Kronig integral relations.
      <p>
      Reference: B. L. Henke, E. M. Gullikson, and J. C. Davis,
      Atomic Data and Nuclear Data Tables Vol. 54 No. 2 (July 1993).<br>
      ftp://grace.lbl.gov/pub/sf/<br>
      See also:
      http://www.esrf.fr/computing/expg/subgroups/theory/DABAX/tmp_file/FileDesc.html
   */
  class Henke {
    public:
      //! Search Henke table for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          See also: eltbx::StripLabel()
       */
      Henke(const std::string& Label, bool Exact = false);
      //! Return scattering factor label.
      inline const char* Label() const { return m_Label_Z_Efpfdp->Label; }
      //! Return atomic number.
      inline int Z() const { return m_Label_Z_Efpfdp->Z; }
      //! Compute f-prime (f') and f-double-prime (f") for given energy.
      /*! f-prime and f-double-prime are determined by linear
          interpolation.<br>
          See also: cctbx::constants::factor_keV_Angstrom
       */
      fpfdp operator()(double Energy);
    private:
      const eltbx::detail::Label_Z_Efpfdp* m_Label_Z_Efpfdp;
  };

} // eltbx

#endif // CCTBX_ELTBX_HENKE_H
