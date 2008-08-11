#ifndef CCTBX_SGTBX_SELECT_GENERATORS_H
#define CCTBX_SGTBX_SELECT_GENERATORS_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace sgtbx { namespace select_generators {

  struct any
  {
    any() : n_gen(0) {}

    any(space_group const& sg, int z2p_r_den, int z2p_t_den);

    std::size_t n_all() const
    {
      if (z_inv_t.is_valid()) return n_gen + 1;
      return n_gen;
    }

    void set_primitive();

    change_of_basis_op z2p_op;
    tr_vec z_inv_t;
    tr_vec p_inv_t;
    std::size_t n_gen;
    rt_mx z_gen[2];
    rt_mx p_gen[2];
  };

  struct standard : any
  {
    standard(space_group const& work_sg,
             int z2p_r_den,
             int z2p_t_den,
             matrix_group::code const& point_group_mx_group_code);
  };

}}} // namespace cctbx::sgtbx::select_generators

#endif // CCTBX_SGTBX_SELECT_GENERATORS_H
