// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ELTBX_CAASF_H
#define CCTBX_ELTBX_CAASF_H

#include <cstddef>
#include <string>
#include <cctbx/fixes/cmath>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/caasf.h>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace eltbx {

  static const float caasf_undefined = -99999.;

  namespace detail {
    template <std::size_t N>
    struct CAASF_Raw {
      const char* Label;
      float a[N], b[N], c;
    };
  }

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Currently used to work with coefficients from the International
      Tables Volume B (1992) (template parameter N = 4) and Waasmaier &
      Kirfel (1995), Acta Cryst. A51, 416-431 (template parameter N = 5).
   */
  template <std::size_t N>
  class CAASF {
    public:
      //! Constructor. For internal use only.
      CAASF(const detail::CAASF_Raw<N>* TableRaw,
            const char* Table,
            const std::string& Label,
            bool Exact)
        : m_Entry(0), m_Table(Table)
      {
        std::string WorkLabel = StripLabel(Label, Exact);
        int m = 0;
        for (const detail::CAASF_Raw<N>* Entry = TableRaw;
             Entry->Label;
             Entry++) {
          int i = eltbx::MatchLabels(WorkLabel, Entry->Label);
          if (i < 0) {
            if (Entry->c == eltbx::caasf_undefined) {
              throw eltbx::error(
                "Analytical approximation for given label not known.");
            }
            m_Entry = Entry;
            return;
          }
          if (i > m && !isdigit(Entry->Label[i - 1])) {
            m = i;
            m_Entry = Entry;
          }
        }
        if (Exact || !m_Entry) {
          throw eltbx::error("Unknown scattering factor label.");
        }
      }
      //! Label of table. Currently either "IT1992" or "WK1995".
      inline const char* Table() const { return m_Table; }
      //! Scattering factor label. E.g. "Si4+".
      inline const char* Label() const { return m_Entry->Label; }
      //! Number of a and b coefficients. Currently either 4 or 5.
      /*! Note that the total number of coefficients is 2*n_ab()+1.
       */
      inline int n_ab() const { return N; }
      //! Return coefficient a(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
       */
      inline float a(int i) const { return m_Entry->a[i]; }
      //! Return coefficient b(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
       */
      inline float b(int i) const { return m_Entry->b[i]; }
      //! Return coefficient c.
      inline float c() const { return m_Entry->c; }
      //! Return the analytical approximation to the scattering factor.
      /*! stol2 = sin-theta-over-lambda-squared: (sin(theta)/lambda)^2<br>
          See also: uctbx::UnitCell::Q()
       */
      double operator()(double stol2)
      {
        double sf = c();
        for (int i = 0; i < n_ab(); i++)
          sf += a(i) * std::exp(-b(i) * stol2);
        return sf;
      }
    private:
      const char *m_Table;
      const detail::CAASF_Raw<N>* m_Entry;
  };

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from the International Tables Volume B (1992).
      <p>
      See also: class CAASF
   */
  class CAASF_IT1992: public CAASF<4> {
    public:
      //! Lookup coefficients for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::StripLabel()
          <p>
          Note that the other methods of this class are inhereted from
          class CAASF.
       */
      CAASF_IT1992(const std::string& Label, bool Exact = false);
  };

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431.
      <p>
      See also: class CAASF
   */
  class CAASF_WK1995: public CAASF<5> {
    public:
      //! Lookup coefficients for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::StripLabel()
          <p>
          Note that the other methods of this class are inhereted from
          class CAASF.
       */
      CAASF_WK1995(const std::string& Label, bool Exact = false);
  };

} // eltbx

#endif // CCTBX_ELTBX_CAASF_H
