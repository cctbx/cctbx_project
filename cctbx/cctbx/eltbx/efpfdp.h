// $Id$

#ifndef CCTBX_ELTBX_EFPFDP_H
#define CCTBX_ELTBX_EFPFDP_H

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
      //! Return f-prime (f').
      inline float fp() const { return m_fp; }
      //! Return f-double-prime (f").
      inline float fdp() const { return m_fdp; }
      //! Test if f-prime is undefined.
      inline bool isValid_fp() const { return m_fp != Efpfdp_undefined; }
      //! Test if f-double-prime is undefined.
      inline bool isValid_fdp() const { return m_fdp != Efpfdp_undefined; }
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
