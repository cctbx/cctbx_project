// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ELTBX_CAASF_H
#define CCTBX_ELTBX_CAASF_H

#include <cstddef>
#include <string>
#include <ctype.h> // cannot use cctype b/o non-conforming compilers
#include <cctbx/fixes/cmath>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/caasf.h>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace cctbx { namespace eltbx {

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
      Tables Volume C (1992) (template parameter N = 4) and Waasmaier &
      Kirfel (1995), Acta Cryst. A51, 416-431 (template parameter N = 5).
   */
  template <std::size_t N>
  class CAASF {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      CAASF() : m_Entry(0), m_Table(0) {}
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
      const char* Table() const { return m_Table; }
      //! Scattering factor label. E.g. "Si4+".
      const char* Label() const { return m_Entry->Label; }
      //! Number of a and b coefficients. Currently either 4 or 5.
      /*! Note that the total number of coefficients is 2*n_ab()+1.
       */
      static std::size_t n_ab() { return N; }
      //! Return coefficient a(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
       */
      float a(int i) const { return m_Entry->a[i]; }
      //! Return coefficient b(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
       */
      float b(int i) const { return m_Entry->b[i]; }
      //! Return coefficient c.
      float c() const { return m_Entry->c; }
      //! Return the analytical approximation to the scattering factor.
      /*! stol2 = sin-theta-over-lambda-squared: (sin(theta)/lambda)^2<br>
          See also: uctbx::UnitCell::Q()
       */
      double stol2(double stol2_value) const
      {
        double sf = c();
        for (int i = 0; i < n_ab(); i++)
          sf += a(i) * std::exp(-b(i) * stol2_value);
        return sf;
      }
      //! Return the analytical approximation to the scattering factor.
      /*! See also: stol2(), uctbx::UnitCell::Q()
       */
      double stol(double stol_value) const {
        return stol2(stol_value * stol_value);
      }
      //! Return the analytical approximation to the scattering factor.
      /*! See also: stol2(), uctbx::UnitCell::Q()
       */
      double Q(double Q_value) const { return stol2(Q_value / 4.); }
    private:
      const char *m_Table;
      const detail::CAASF_Raw<N>* m_Entry;
  };

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from the International Tables for Crystallography,
      Volume C, Mathematical, Physical and Chemical Tables,
      Edited by A.J.C. Wilson, Kluwer Academic Publishers,
      Dordrecht/Boston/London, 1992.
      <p>
      Table 6.1.1.4 (pp. 500-502):
      Coefficients for analytical approximation to the scattering factors
      of Tables 6.1.1.1 and 6.1.1.3.
      <p>
      Table 6.1.1.4 is a reprint of Table 2.2B, pp. 99-101,
      International Tables for X-ray Crystallography, Volume IV,
      The Kynoch Press: Birmingham, England, 1974.
      There is just one difference for <tt>Tl3+</tt>:<pre>
                IT Vol IV 1974: a2 = 18.3841
                IT Vol C  1992: a2 = 18.3481</pre>
      <p>
      The value from IT Vol IV is presumed to be the correct one
      and used here.
      <p>
      If possible, the class CAASF_WK1995 should be used
      instead of this class.
      The coefficients of Waasmaier & Kirfel give more precise
      approximations than the coefficients of Volume C
      and some errors are corrected.
      <p>
      See also: class CAASF, class CAASF_WK1995
   */
  class CAASF_IT1992: public CAASF<4> {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      CAASF_IT1992() {}
      //! Lookup coefficients for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::StripLabel()
          <p>
          Note that the other methods of this class are inherited from
          class CAASF.
       */
      CAASF_IT1992(const std::string& Label, bool Exact = false);
  };

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431.
      <p>
      The original data can be found at:
      <pre>
      ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat
      </pre>
      File picked up Jul  4, 1995.
      File verified  Sep 12, 2001.
      <p>
      Header of "sfac.dat":
      <pre>
        - Fitparameters of all atoms/ions (with the excepetion of O1-)
          from publication "New Analytical Scattering Factor Functions
          for Free Atoms and Ions", D. Waasmaier & A. Kirfel, Acta Cryst.
          A 1995, in press)

        - Fit for O1- based on the tabulated values of Table 2 (D.REZ,
          P.REZ & I.GRANT, Acta Cryst. (1994), A50, 481-497).

        - Fits for all other atoms/ions based on the tabulated values
          of Table 6.1.1.1 (atoms) and Table 6.1.1.3 (ions)
          (International Tables for Crystallography, Vol. C, 1992).
      </pre>
      <p>
      See also: class CAASF, class CAASF_IT1992
   */
  class CAASF_WK1995: public CAASF<5> {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      CAASF_WK1995() {}
      //! Lookup coefficients for the given scattering factor label.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::StripLabel()
          <p>
          Note that the other methods of this class are inherited from
          class CAASF.
       */
      CAASF_WK1995(const std::string& Label, bool Exact = false);
  };

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_CAASF_H
