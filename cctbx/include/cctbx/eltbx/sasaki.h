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

namespace cctbx { namespace eltbx {

  namespace detail { namespace sasaki {

    const long n_raw = 280; // All tables have exactly 280 data points.
    const double first_wide = 0.1; // All tables in fpwide.tbl start at 0.1.
    const double wide_incr = 0.01;
    const double edge_incr = 0.0001;

    struct raw {
      float fp;
      float fdp;
    };

    struct info {
      char* Label;
      int Z;
      raw* wide;
      double first_k;
      raw* k;
      double first_l1;
      raw* l1;
      double first_l2;
      raw* l2;
      double first_l3;
      raw* l3;
    };

  }} // namespace detail::sasaki

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
      //! Default constructor. Calling certain methods may cause crashes!
      Sasaki() : m_info(0) {}
      //! Search Sasaki table for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          See also: eltbx::StripLabel()
       */
      explicit
      Sasaki(const std::string& Label, bool Exact = false);
      //! Return scattering factor label.
      const char* Label() const { return m_info->Label; }
      //! Return atomic number.
      int Z() const { return m_info->Z; }
      //! Compute f-prime (f') and f-double-prime (f") for given energy [eV].
      /*! f-prime and f-double-prime are determined by linear
          interpolation.<br>
          See also: cctbx::constants::factor_keV_Angstrom
       */
      fpfdp operator()(double Energy);
    private:
      const detail::sasaki::info* m_info;
  };

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_SASAKI_H
