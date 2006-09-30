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
