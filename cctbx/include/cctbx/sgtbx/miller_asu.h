// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 27: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#ifndef CCTBX_SGTBX_MILLER_ASU_H
#define CCTBX_SGTBX_MILLER_ASU_H

#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_ops.h>
#include <cctbx/sgtbx/miller_ref_asu.h>

namespace cctbx {
  namespace sgtbx {

    //! Access to general contiguous reciprocal space asymmetric units.
    /*! class ReferenceReciprocalSpaceASU implements 12 contiguous
        reciprocal space asymmetric units that cover the 230
        reference settings. The algorithm for the determination of the
        space group type (class SpaceGroupInfo) is used to derive a
        change-of-basis matrix for the transformation of the tabulated
        asymmetric units. In this way a contiguous asymmetric unit is
        available for any arbitrary setting.
     */
    class ReciprocalSpaceASU
    {
      public:
        //! Default constructor.
        /*! Default-constructed instances will throw exceptions if
            some of the member functions are used.
         */
        ReciprocalSpaceASU()
          : m_CBOp(), m_isReferenceASU(true), m_ReferenceASU() {}
        //! Initialization.
        /*! Based on the space group number (SpaceGroupInfo::SgNumber()),
            the Laue group class is derived which is in turn used
            to select the corresponding tabulated
            ReferenceReciprocalSpaceASU.
         */
        ReciprocalSpaceASU(const SpaceGroupInfo& SgInfo);
        //! Access to the selected tabulated ReferenceReciprocalSpaceASU.
        const ReferenceReciprocalSpaceASU* ReferenceASU() const {
          return m_ReferenceASU;
        }
        //! Access to the change-of-basis operator.
        /*! This operator is a copy of SgInfo.CBOp() as passed to
            the constructor.
         */
        const ChOfBasisOp& CBOp() const { return m_CBOp; }
        //! Test if the given asymmetric unit is one of the tabulated units.
        /*! This test is equivalent to the test CBOp().M().Rpart().isUnit().
            That is, it is tested if the rotation part of the
            change-of-basis matrix is the unit matrix. (If this
            is the case, some optimizations are possible.)
         */
        bool isReferenceASU() const { return m_isReferenceASU; }
        //! Test if the given Miller index is in the asymmetric unit.
        /*! The change-of-basis matrix is used to transform
            the Miller index (H * CBOp().InvM()). It is then
            tested if the result is in the tabulated reference
            asymmetric unit.
         */
        bool isInASU(const miller::Index& H) const {
          if (m_isReferenceASU) return m_ReferenceASU->isInASU(H);
          return m_ReferenceASU->isInASU(H * m_CBOp.InvM().Rpart());
        }
        int asu_sign(const miller::Index& h,
                     const miller::Index& minus_h) const {
          if (m_isReferenceASU) return m_ReferenceASU->asu_sign(h, minus_h);
          miller::Index h_ref = h * m_CBOp.InvM().Rpart();
          return m_ReferenceASU->asu_sign(h_ref);
        }
        int asu_sign(const miller::Index& h) const {
          if (m_isReferenceASU) return m_ReferenceASU->asu_sign(h);
          miller::Index h_ref = h * m_CBOp.InvM().Rpart();
          return m_ReferenceASU->asu_sign(h_ref);
        }
      private:
        ChOfBasisOp m_CBOp;
        bool m_isReferenceASU;
        const ReferenceReciprocalSpaceASU* m_ReferenceASU;
    };

  } // namespace sgtbx
} // namespace cctbx

#endif // CCTBX_SGTBX_MILLER_ASU_H
