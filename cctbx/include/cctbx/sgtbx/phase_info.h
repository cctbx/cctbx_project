// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Renamed from miller.h to phase_info.h (rwgk)
     2001 Oct 27: Redesign: AsymIndex (rwgk)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_PHASE_INFO_H
#define CCTBX_SGTBX_PHASE_INFO_H

#include <cctbx/miller.h>

namespace cctbx { namespace sgtbx {

  class SpaceGroup; // forward declaration

  /*! \brief Handling of phase restrictions and optional evaluation
      of conditions for systematically absent reflections.
   */
  /*! A reflection with the Miller index H is "centric" if
      there is a symmetry operation with rotation part R such that
      H*R = -H. The phase of a centric reflection is restricted to
      two phase angels (modulo pi).
      <p>
      A reflection with the Miller index H is "systematically absent"
      if there is a symmetry operation with the rotation part R and
      the translation part T such that H*R == H and H*T != 0 mod 1.
   */
  class PhaseInfo
  {
    public:
      //! Default constructor. Some data members are not initialized!
      PhaseInfo() {}

      //! Determination of the phase restriction for a given Miller index.
      /*! If no_test_sys_absent == false, it is also tested if a
          reflection with the given Miller index is systematically
          absent. If no_test_sys_absent == true, a faster algorithm
          is used that only determines the phase restriction. In the
          latter case the isSysAbsent() member function must not
          be used (an exception is thrown otherwise).
       */
      PhaseInfo(SpaceGroup const& sgops, miller::Index const& h,
                bool no_test_sys_absent = false);

      //! Initialization with known product H*T and given base factor.
      /*! sys_abs_was_tested indicates if it is known if a reflection
          with the given Miller index is systematically absent.
          The isSysAbsent() member function must not be used
          if sys_abs_was_tested == false (an exception is thrown
          otherwise).
       */
      PhaseInfo(int HT, int TBF, bool sys_abs_was_tested)
        : m_HT(HT), m_TBF(TBF), m_SysAbsWasTested(sys_abs_was_tested)
      {}

      //! Test if isSysAbsent() can be used.
      bool SysAbsWasTested() const { return m_SysAbsWasTested; }

      //! Test for systematically absent reflection.
      /*! Use SysAbsWasTested() to determine if this test can be used.
          If the result of SysAbsWasTested() == false, an exception
          is thrown if this test is used.
       */
      bool isSysAbsent() const
      {
        cctbx_assert(m_SysAbsWasTested);
        return m_HT == -2;
      }

      //! Test if there actually is a phase restriction.
      /*! See class details.
       */
      bool isCentric() const { return m_HT >= 0; }

      //! Phase shift H*T (mod 1) corresponding to H*R = -H.
      /*! Low-level information for computing the restricted phases.
          HT() is multiplied by a base factor TBF() in order to obtain
          an integer value.
          <p>
          See also: HT_angle()
       */
      int HT() const { return m_HT; }
      //! Translation base factor.
      /*! This is the factor by which HT() is multiplied.
       */
      int TBF() const { return m_TBF; }

      //! Phase restriction in radians or degrees.
      /*! The return value is -1 if the phase is not restricted,
          and >= 0 and < pi or 180 otherwise.
       */
      double HT_angle(bool deg = false) const
      {
        if (!isCentric()) return -1.;
        return (ht_period(deg) * m_HT) / m_TBF;
      }

      /*! \brief Test if phase phi (with given Period) is
          compatible with restriction.
       */
      //! Test if phase phi is compatible with restriction.
      /*! The tolerance compensates for rounding errors.
       */
      bool isValidPhase(
        double phi, bool deg = false, double tolerance = 1.e-5) const;

      //! Nearest valid phase.
      /*! For acentric reflections equivalent to the input phase phi.
          For centric reflections, the restricted phase which is
          closest to the input phase phi.
       */
      double nearest_valid_phase(double phi, bool deg = false) const;

    private:
      double ht_period(bool deg) const
      {
        if (deg) return 180.;
        return cctbx::constants::pi;
      }

      int m_HT;
      int m_TBF;
      bool m_SysAbsWasTested;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_PHASE_INFO_H
