/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored (rwgk)
     2001 Sep: created from fragments in type.cpp, seminvariant.cpp (rwgk)
 */

#include <cctbx/sgtbx/select_generators.h>
#include <cctbx/sgtbx/rot_mx_info.h>

namespace cctbx { namespace sgtbx { namespace select_generators {

  any::any(space_group const& sg,
           int z2p_r_den,
           int z2p_t_den)
  : n_gen(0)
  {
    using namespace crystal_system;
    using scitbx::fn::absolute;

    z2p_op = sg.z2p_op(z2p_r_den, z2p_t_den);

    z_inv_t = sg.inv_t(true);
    p_inv_t = tr_vec(0);
    for(std::size_t i=0;i<2;i++) z_gen[i] = rt_mx(0, 0);
    for(std::size_t i=0;i<2;i++) p_gen[i] = rt_mx(0, 0);

    int principal_proper_order = 0;

    matrix_group::code point_group_mx_group_code = sg.point_group_type();
    switch (point_group_mx_group_code.crystal_system())
    {
      case triclinic:
        break;

      case monoclinic:
        z_gen[0] = sg.smx(1);
        n_gen = 1;
        break;

      case orthorhombic:
        z_gen[0] = sg.smx(1);
        z_gen[1] = sg.smx(2);
        n_gen = 2;
        break;

      case tetragonal:
                                     principal_proper_order = 4;
      case trigonal:
        if (!principal_proper_order) principal_proper_order = 3;
      case hexagonal:
        if (!principal_proper_order) principal_proper_order = 6;
        {
          rot_mx_info principal_ri;
          std::size_t i=1;
          for(;i<sg.n_smx();i++) {
            principal_ri = sg.smx(i).r().info();
            if (absolute(principal_ri.type()) == principal_proper_order) {
              if (principal_ri.sense() > 0) {
                z_gen[0] = sg.smx(i);
                n_gen++;
                break;
              }
            }
          }
          CCTBX_ASSERT(n_gen == 1);
          std::size_t i_principal = i;
          for(i=1;i<sg.n_smx();i++) {
            if (i == i_principal) continue;
            rot_mx_info ri(sg.smx(i).r());
            if (absolute(ri.type()) == 2) {
              if (principal_ri.ev() != ri.ev()) {
                z_gen[1] = sg.smx(i);
                n_gen++;
                break;
              }
            }
          }
        }
        break;

      case cubic:
        for(std::size_t i=1;i<sg.n_smx();i++) {
          rot_mx_info ri(sg.smx(i).r());
          if      (absolute(ri.type()) == 3) {
            if (ri.sense() > 0) {
              if (!z_gen[0].is_valid()) {
                z_gen[0] = sg.smx(i);
                n_gen++;
                if (n_gen == 2) break;
              }
            }
          }
          else if (absolute(ri.type()) == sg.n_smx() / 6) {
            if (ri.sense() >= 0) {
              if (!z_gen[1].is_valid()) {
                z_gen[1] = sg.smx(i);
                n_gen++;
                if (n_gen == 2) break;
              }
            }
          }
        }
        CCTBX_ASSERT(n_gen == 2);
        break;

      default:
        throw CCTBX_INTERNAL_ERROR();
    }
  }

  void any::set_primitive()
  {
    for (std::size_t i=0;i<n_gen;i++) {
      p_gen[i] = z2p_op(z_gen[i]).mod_positive();
    }
    if (z_inv_t.is_valid()) {
      p_inv_t = z2p_op(z_inv_t, -1).mod_positive();
    }
  }

  standard::standard(space_group const& work_sg,
                     int z2p_r_den,
                     int z2p_t_den,
                     matrix_group::code const& point_group_mx_group_code)
  {
    using namespace crystal_system;
    using scitbx::fn::absolute;

    const sg_vec3 ev_001( 0, 0, 1);
    const sg_vec3 ev_100( 1, 0, 0);
    const sg_vec3 ev_110( 1, 1, 0);
    const sg_vec3 ev_m10(-1, 1, 0);
    const sg_vec3 ev_111( 1, 1, 1);

    z2p_op = work_sg.z2p_op(z2p_r_den, z2p_t_den);

    z_inv_t = work_sg.inv_t(true);
    p_inv_t = tr_vec(0);
    for(std::size_t i=0;i<2;i++) z_gen[i] = rt_mx(0, 0);
    for(std::size_t i=0;i<2;i++) p_gen[i] = rt_mx(0, 0);

    int principal_proper_order = 0;

    switch (point_group_mx_group_code.crystal_system())
    {
      case triclinic:
        break;

      case monoclinic:
        z_gen[0] = work_sg.smx(1);
        n_gen = 1;
        break;

      case orthorhombic:
        for(std::size_t i=1;i<work_sg.n_smx();i++) {
          rot_mx_info ri(work_sg.smx(i).r());
          if      (ri.ev() == ev_001) {
            z_gen[0] = work_sg.smx(i); n_gen++;
          }
          else if (ri.ev() == ev_100) {
            z_gen[1] = work_sg.smx(i); n_gen++;
          }
        }
        CCTBX_ASSERT(n_gen == 2);
        break;

      case tetragonal:
                                     principal_proper_order = 4;
      case trigonal:
        if (!principal_proper_order) principal_proper_order = 3;
      case hexagonal:
        if (!principal_proper_order) principal_proper_order = 6;

        for(std::size_t i=1;i<work_sg.n_smx();i++) {
          rot_mx_info ri(work_sg.smx(i).r());
          if (absolute(ri.type()) == principal_proper_order) {
            if (ri.sense() > 0) {
              z_gen[0] = work_sg.smx(i); n_gen++;
            }
          }
          else if (principal_proper_order == 4) {
            if (ri.ev() == ev_100) {
              z_gen[1] = work_sg.smx(i); n_gen++;
            }
          }
          else if (principal_proper_order == 3) {
            if      (ri.ev() == ev_m10) {
              z_gen[1] = work_sg.smx(i); n_gen++;
            }
            else if (ri.ev() == ev_110) {
              z_gen[1] = work_sg.smx(i); n_gen++;
            }
          }
          else { // principal_proper_order == 6
            if (ri.ev() == ev_m10) {
              z_gen[1] = work_sg.smx(i); n_gen++;
            }
          }
        }
        CCTBX_ASSERT(n_gen == 1 || n_gen == 2);
        for (std::size_t i=0;i<n_gen;i++) CCTBX_ASSERT(z_gen[i].is_valid());
        break;

      case cubic:
        for(std::size_t i=1;i<work_sg.n_smx();i++) {
          rot_mx_info ri(work_sg.smx(i).r());
          if      (absolute(ri.type()) == 4) {
            if (ri.sense() > 0 && ri.ev() == ev_001) {
              if (!z_gen[0].is_valid()) n_gen++;
              z_gen[0] = work_sg.smx(i);
            }
          }
          else if (absolute(ri.type()) == 2) {
            if (!z_gen[0].is_valid() && ri.ev() == ev_001) {
              z_gen[0] = work_sg.smx(i); n_gen++;
            }
          }
          else if (absolute(ri.type()) == 3) {
            if (ri.sense() > 0 && ri.ev() == ev_111) {
              CCTBX_ASSERT(!z_gen[1].is_valid());
              z_gen[1] = work_sg.smx(i); n_gen++;
            }
          }
        }
        CCTBX_ASSERT(n_gen == 1 || n_gen == 2);
        for (std::size_t i=0;i<n_gen;i++) CCTBX_ASSERT(z_gen[i].is_valid());
        break;

      default:
        throw CCTBX_INTERNAL_ERROR();
    }

    // Tidy generators
    if (z_inv_t.is_valid()) {
      for (std::size_t i=0;i<n_gen;i++) {
        if (z_gen[i].r().num().determinant() < 0) {
          z_gen[i] = z_gen[i].pre_multiply_inv_t(z_inv_t);
        }
      }
    }
    for (std::size_t i=0;i<n_gen;i++) z_gen[i].mod_positive_in_place();
  }

}}} // namespace cctbx::sgtbx::select_generators
