// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: split of miller_asu.h ((R.W. Grosse-Kunstleve)
     2001 Oct 27: Redesign: AsymIndex (rwgk)
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#ifndef CCTBX_SGTBX_MILLER_REF_ASU_H
#define CCTBX_SGTBX_MILLER_REF_ASU_H

#include <cctbx/miller.h>
#include <cctbx/sgtbx/group_codes.h>

namespace cctbx { namespace sgtbx {

  /*! \brief Contiguous reciprocal space asymmetric units for
      the 230 reference settings.
   */
  /*! 12 contiguous reciprocal space asymmetric units (11 Laue
      classes, two settings for Laue class -3m) are
      tabulated. The tabulated asymmetric units are
      compatible with the asymmetric units of the CCP4
      package.
      <p>
      This implementation is based on work by
      <a href="http://www.ysbl.york.ac.uk/~cowtan/"
      >Kevin Cowtan</a>.
      <p>
      See also: class ReciprocalSpaceASU, class MillerIndexGenerator
   */
  class ReferenceReciprocalSpaceASU {
    public:
      //! Returns one of exactly 12 Laue group codes.
      /*! The labels of the possible return codes are:<br>
            -1, 2/m, mmm, 4/m, 4/mmm, -3, -31m, -3m1, 6/m, 6/mmm, m-3, m-3m
          <p>
          For Laue class -3m there are two possible orientations of the
          mirror plane with respect to the periodic lattice.
       */
      virtual tables::MatrixGroup::Code LaueGroupCode() const {
        throw cctbx_internal_error();
      }
      //! Test if given Miller index is in the tabulated asymmetric unit.
      virtual bool isInASU(const Miller::Index& H) const {
        throw cctbx_internal_error();
      }
      //! XXX
      int asu_sign(const Miller::Index& h,
                   const Miller::Index& minus_h) const {
        if      (isInASU(      h)) return  1;
        else if (isInASU(minus_h)) return -1;
        return 0;
      }
      //! XXX
      int asu_sign(const Miller::Index& h) const {
        return asu_sign(h, -h);
      }
      //! String representation of the tabluated asymmetric unit.
      /*! Example: "h>=k and k>=0 and (k>0 or l>=0)"
       */
      virtual const char* representation() const {
        throw cctbx_internal_error();
      }
      //! "Cut parameters" for building Miller indices.
      /*! When building (or generating) a large list of Miller indices,
          it is useful to restrict the loop over all possible indices
          to 1/2, 1/4, or 1/8 of reciprocal space, if possible.
          <p>
          The cut parameters are used in the next() method
          of the class MillerIndexGenerator. In general it
          should be much more convenient to use that higher-level
          class rather than hand-crafting a loop for building
          Miller indices.
          <p>
          Each element of CutP is either -1 or 0. A value
          of 0 indicates that the corresponding negative half-space
          can be omitted in the loop over possible indices.
          <p>
          Friedel symmetry is implied. If the Friedel mates
          are needed explicitly, they have to be added in a
          separate step. Note that the Friedel mate appears
          explicitly only for acentric reflections (use e.g.
          !SpaceGroup::isCentric(H) to determine which reflections
          are acentric).
       */
      virtual const af::int3& getCutParameters() const {
        throw cctbx_internal_error();
      }
  };

  bool isInReferenceReciprocalSpaceASU_1b(Miller::Index const& h);

  const ReferenceReciprocalSpaceASU*
  LookupReferenceReciprocalSpaceASU(tables::MatrixGroup::Code group_code);

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_MILLER_REF_ASU_H
