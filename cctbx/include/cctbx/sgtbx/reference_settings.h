/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/sgtbx/reference.h (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: Created, based on sglite/sgrefset.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_REFERENCE_SETTINGS_H
#define CCTBX_SGTBX_REFERENCE_SETTINGS_H

#include <cctbx/sgtbx/group_codes.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace sgtbx { namespace reference_settings {

  /* Hall symbols for the reference settings of the 230 crystallographic
     space groups.

     The reference settings chosen are identical to those listed in
     International Tables for Crystallography Vol. A. For the cases
     where more than one setting is given in the International Tables,
     the following choices have been made:
       - For monoclinic space groups: unique axis b and cell choice 1.
       - For space groups with two origin choices: origin choice 2.
       - Rhombohedral space groups: hexagonal axes.
   */
  const char*
  hall_symbol_table(std::size_t i);

  /* Matrix group codes (Boisen & Gibbs, 1990, pp. 225-228)
     corresponding to the reference settings above.
   */
  matrix_group::code
  const& matrix_group_code_table(std::size_t i);

  namespace normalizer {

    struct addl_generators
    {
      const char* k2l; // operations which generate L from K
      const char* l2n; // operations which generate N from L
    };

    /* Table of 'additional generators' of the Euclidean and affine
       normalizers, corresponding to the reference settings above.
       Reference: Int. Tab. Vol. A Section 15.3
     */
    addl_generators const&
    addl_generators_table(std::size_t i);

    af::shared<rt_mx>
    get_addl_generators(
      int sg_number,
      bool flag_affine,
      bool flag_k2l,
      bool flag_l2n);

    bool
    check_monoclinic_affine_restrictions(
      int sg_number,
      rot_mx const& r);

    void
    get_monoclinic_affine_trial_ranges(
      rot_mx const& cb_mx_r,
      int& r00,
      int& r22);

  } // namespace normalizer

  namespace wyckoff {

    struct raw_position
    {
      int m;
      const char *xyz;
    };

    struct raw_table
    {
      int n;
      const raw_position *op;
    };

    const raw_table&
    raw_tables(std::size_t i);

    int
    general_position_multiplicities(std::size_t i);

  } // namespace wyckoff

}}} // namespace cctbx::sgtbx::reference_settings

#endif // CCTBX_SGTBX_REFERENCE_SETTINGS_H
