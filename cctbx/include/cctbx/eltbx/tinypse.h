// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ELTBX_TINYPSE_H
#define CCTBX_ELTBX_TINYPSE_H

#include <string>

namespace cctbx { namespace eltbx {

  namespace detail {
    struct TinyPSE_RawEntry {
      int Z;
      const char* Symbol;
      const char* Name;
      float Weight;
    };
  }

  //! Tiny Periodic System of Elements.
  /*! Reference:<br>
      CRC Handbook of Chemistry & Physics, 63rd edition, 1982-1983<br>
      CRC Handbook of Chemistry & Physics, 70th edition, 1989-1990
   */
  class TinyPSE {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      TinyPSE() : m_RawEntry(0) {}
      //! Lookup table entry by element symbol.
      /*! If Exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          See also: eltbx::StripLabel()
       */
      explicit
      TinyPSE(const std::string& Label, bool Exact = false);
      //! Lookup table entry by atomic number.
      TinyPSE(int Z);
      //! Return atomic number.
      int Z() const { return m_RawEntry->Z; }
      //! Return element symbol.
      const char* Symbol() const { return m_RawEntry->Symbol; }
      //! Return element name.
      const char* Name() const { return m_RawEntry->Name; }
      //! Return atomic weight.
      float Weight() const { return m_RawEntry->Weight; }
    private:
      const detail::TinyPSE_RawEntry* m_RawEntry;
  };

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_TINYPSE_H
