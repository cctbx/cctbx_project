// $Id$

#ifndef CCTBX_ELTBX_HENKE_H
#define CCTBX_ELTBX_HENKE_H

#include <string>
#include <cctbx/eltbx/efpfdp.h>

namespace eltbx {

  //! Access to Henke tables.
  /*! Reference: B. L. Henke, E. M. Gullikson, and J. C. Davis,
      Atomic Data and Nuclear Data Tables Vol. 54 No. 2 (July 1993).<br>
      ftp://grace.lbl.gov/pub/sf/
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
