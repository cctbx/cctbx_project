// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#ifndef CCTBX_ELTBX_SASAKI_H
#define CCTBX_ELTBX_SASAKI_H

#include <string>
#include <cctbx/eltbx/efpfdp.h>

namespace eltbx {

  //! Access to Sasaki tables.
  /*! Sasaki tables are available for elements with Z=4-83 and Z=92.
      They are valid in the energy range 4-124 keV and they have a fine
      step size close to the absorption edges (K,L1,L2,L3). The tables
      are therefore suitable for use in connection with anomalous
      diffraction experiments.
      <p>
      Reference: S.Sasaki (1989) Numerical Tables of Anomalous
      Scattering Factors Calculated by the Cromer and Liberman Method,
      KEK Report, 88-14, 1-136<br>
      ftp://pfweis.kek.jp/pub/Sasaki-table/<br>
      See also: http://www.esrf.fr/computing/expg/subgroups/theory/DABAX/tmp_file/FileDesc.html
   */
  class Sasaki {
    public:
      //! Search Sasaki table for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          See also: eltbx::StripLabel()
       */
      Sasaki(const std::string& Label, bool Exact = false);
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

#endif // CCTBX_ELTBX_SASAKI_H
