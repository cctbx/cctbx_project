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

namespace eltbx {

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

  static const float Efpfdp_undefined = -9999.00;

  //! Helper class for passing f' (f-prime) and f" (f-double-prime).
  class fpfdp {
    public:
      //! Constructor. For internal use only.
      inline fpfdp(float fp, float fdp) : m_fp(fp), m_fdp(fdp) {}
      //! Test if f-prime is defined.
      inline bool isValid_fp() const { return m_fp != Efpfdp_undefined; }
      //! Test if f-double-prime is defined.
      inline bool isValid_fdp() const { return m_fdp != Efpfdp_undefined; }
      //! Test if f-prime and f-double-prime are defined.
      inline bool isValid() const { return isValid_fp() && isValid_fdp(); }
      //! Return f-prime (f').
      inline float fp() const { return m_fp; }
      //! Return f-double-prime (f").
      inline float fdp() const { return m_fdp; }
      //! Return complex(f-prime, f-double-prime).
      inline std::complex<float> operator()() const {
        return std::complex<float>(m_fp, m_fdp);
      }
    private:
      float m_fp;
      float m_fdp;
  };

  namespace detail {
    const Label_Z_Efpfdp* FindEntry(const Label_Z_Efpfdp* Tables,
                                    const std::string& WorkLabel,
                                    bool Exact);
    fpfdp interpolate(const Label_Z_Efpfdp* m_Label_Z_Efpfdp, double Energy);
  }

} // eltbx

#endif // CCTBX_ELTBX_EFPFDP_H
