// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 31: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/sgtbx/reference.h>

namespace cctbx {
  namespace sgtbx {

    ReciprocalSpaceASU::ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo)
      : m_CBOp(SgInfo.CBOp()),
        m_isReferenceASU(SgInfo.CBOp().M().Rpart().isUnit()),
        m_ReferenceASU(LookupReferenceReciprocalSpaceASU(
          tables::ReferenceSettings::MatrixGroupCodes[SgInfo.SgNumber()]))
    {}


  } // namespace sgtbx
} // namespace cctbx
