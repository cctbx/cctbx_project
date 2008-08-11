#ifndef CCTBX_SGTBX_REFERENCE_SETTINGS_H
#define CCTBX_SGTBX_REFERENCE_SETTINGS_H

#include <cctbx/sgtbx/group_codes.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace sgtbx {

/*! \brief Tables for the "reference settings" of the
    230 crystallographic space groups.
 */
/*! The reference settings chosen are identical to those listed in
    International Tables for Crystallography Vol. A. For the cases
    where more than one setting is given in the International Tables,
    the following choices have been made:
      - For monoclinic space groups: unique axis b and cell choice 1.
      - For space groups with two origin choices: origin choice 2.
      - Rhombohedral space groups: hexagonal axes.
 */
namespace reference_settings {

  //! Hall symbols for the reference settings.
  const char*
  hall_symbol_table(std::size_t i);

  //! Hermann-Mauguin symbols for the reference settings.
  const char*
  hermann_mauguin_symbol_table(std::size_t i);

  /*! \brief Matrix group codes (Boisen & Gibbs, 1990, pp. 225-228)
      corresponding to the reference settings.
   */
  matrix_group::code
  const& matrix_group_code_table(std::size_t i);

  //! Tables for the Euclidean and affine normalizers.
  namespace normalizer {

    struct addl_generators
    {
      const char* k2l; // operations which generate L from K
      const char* l2n; // operations which generate N from L
    };

    /*! \brief Table of "additional generators" of the Euclidean and affine
        normalizers for to the reference settings.
     */
    /*! Reference: Int. Tab. Vol. A Section 15.3
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

  //! Wyckoff tables.
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
