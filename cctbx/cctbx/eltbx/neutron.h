// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#ifndef CCTBX_ELTBX_NEUTRON_H
#define CCTBX_ELTBX_NEUTRON_H

#include <string>
#include <complex>

namespace cctbx { namespace eltbx {

  namespace detail {
    struct RawNeutronNews1992Record {
      const char* Symbol;
      float       BoundCohScattLengthReal;
      float       BoundCohScattLengthImag;
      float       AbsCrossSect; // For 2200 m/s neutrons
    };
  }

  //! Access to neutron bound scattering lengths & cross-sections.
  /*! Reference:<br>
      Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
      <p>
      http://www.ncnr.nist.gov/resources/n-lengths/list.html
   */
  class NeutronNews1992Record {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      NeutronNews1992Record() : m_RawEntry(0) {}
      //! Search internal table for the given element label.
      /*! If Exact == true, the element label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          See also: eltbx::StripLabel()
       */
      explicit
      NeutronNews1992Record(const std::string& Label, bool Exact = false);
      //! Return element label from internal table.
      const char* Symbol() const {
        return m_RawEntry->Symbol; }
      //! Return bound coherent scattering length (fm) as a complex number.
      /*! 1 fm = 1e-15 m
       */
      std::complex<float> BoundCohScattLength() const {
        return std::complex<float>(m_RawEntry->BoundCohScattLengthReal,
                                   m_RawEntry->BoundCohScattLengthImag);
      }
      //! Return real part of bound coherent scattering length (fm).
      /*! 1 fm = 1e-15 m
       */
      float BoundCohScattLengthReal() const {
        return m_RawEntry->BoundCohScattLengthReal; }
      //! Return imaginary part of bound coherent scattering length (fm).
      /*! 1 fm = 1e-15 m
       */
      float BoundCohScattLengthImag() const {
        return m_RawEntry->BoundCohScattLengthImag; }
      //! Return absorption cross section (barn) for 2200 m/s neutrons.
      /*! 1 barn = 1e-24 cm^2
       */
      float AbsCrossSect() const {
        return m_RawEntry->AbsCrossSect; }
    private:
      const detail::RawNeutronNews1992Record* m_RawEntry;
  };

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_NEUTRON_H
