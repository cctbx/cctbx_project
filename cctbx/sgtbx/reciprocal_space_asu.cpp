/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Renamed miller_asu.cpp -> reciprocal_space_asu.cpp (rwgk)
     2001 Oct: Redesign: AsymIndex (rwgk)
     2001 Sep: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Aug: Redesign of Kevin Cowtan's implementation for the
               handling of CCP4 reciprocal-space asymmetric units.
               Motivation: implementation of MillerIndexGenerator (rwgk).
 */

#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/reference_settings.h>

namespace cctbx { namespace sgtbx { namespace reciprocal_space {

  asu::asu(space_group_type const& sg_type)
  : cb_op_(sg_type.cb_op()),
    is_reference_(sg_type.cb_op().is_identity_op()),
    reference_(lookup_reference_asu(
      reference_settings::matrix_group_code_table(sg_type.number())))
  {}

}}} // namespace cctbx::sgtbx::reciprocal_space
