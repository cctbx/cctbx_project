// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#ifndef CCTBX_ELTBX_EFPFDP_H
#define CCTBX_ELTBX_EFPFDP_H

#include <complex>
#include <string>
#include <ctype.h> // cannot use cctype b/o non-conforming compilers

namespace cctbx { namespace eltbx {

  namespace detail {

    struct Efpfdp {
      float E;
      float fp;
      float fdp;
    };

    struct Label_Z_Efpfdp {
      const char* Label;
      int Z;
      const Efpfdp* Data;
    };

  } // namespace detail

  //! Indicator for undefined values.
  static const float Efpfdp_undefined = -9999.00;

  //! Helper class for passing f' (f-prime) and f" (f-double-prime).
  class fpfdp {
    public:
      //! Constructor.
      explicit
      fpfdp(float fp  = Efpfdp_undefined,
            float fdp = Efpfdp_undefined) : m_fp(fp), m_fdp(fdp) {}
      //! Test if f-prime is defined.
      bool isValid_fp() const { return m_fp != Efpfdp_undefined; }
      //! Test if f-double-prime is defined.
      bool isValid_fdp() const { return m_fdp != Efpfdp_undefined; }
      //! Test if f-prime and f-double-prime are defined.
      bool isValid() const { return isValid_fp() && isValid_fdp(); }
      //! Return f-prime (f').
      float fp() const { return m_fp; }
      //! Return f-double-prime (f").
      float fdp() const { return m_fdp; }
      //! Return complex(f-prime, f-double-prime).
      std::complex<float> operator()() const {
        return std::complex<float>(m_fp, m_fdp);
      }
    private:
      float m_fp;
      float m_fdp;
  };

  namespace detail {

    template <typename TableType>
    const TableType*
    FindEntry(const TableType* Tables,
              const std::string& WorkLabel,
              bool Exact)
    {
      // map D to H
      std::string WL = WorkLabel;
      if (WL == "D") WL = "H";
      int m = 0;
      const TableType* mEntry = 0;
      for (const TableType* Entry = Tables; Entry->Label; Entry++)
      {
        int i = MatchLabels(WL, Entry->Label);
        if (i < 0) return Entry;
        if (i > m && !isdigit(Entry->Label[i - 1])) {
          m = i;
          mEntry = Entry;
        }
      }
      if (Exact || !mEntry) {
        throw error("Unknown scattering factor label.");
      }
      return mEntry;
    }

    fpfdp interpolate(const Label_Z_Efpfdp* m_Label_Z_Efpfdp, double Energy);

  } // namespace detail

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_EFPFDP_H
