#include <cctbx/crystal/symmetry.h>
#include <cctbx/uctbx/fast_minimum_reduction.h>
#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/sgtbx/find_affine.h>
#include <cctbx/sgtbx/rot_mx_info.h>
#include <scitbx/array_family/sort.h>
#include <cstdio>

namespace cctbx { namespace {

  // cctbx.crystal.symmetry.change_of_basis_op_to_minimum_cell()
  sgtbx::change_of_basis_op
  get_cb_op_inp_minimum(
    crystal::symmetry const& input_symmetry)
  {
    sgtbx::change_of_basis_op z2p_op = input_symmetry.space_group().z2p_op();
    uctbx::unit_cell p_cell = input_symmetry.unit_cell().change_basis(
      z2p_op.c_inv().r().as_double());
    uctbx::fast_minimum_reduction<> red(p_cell);
    return sgtbx::change_of_basis_op(
      sgtbx::rt_mx(sgtbx::rot_mx(red.r_inv(), 1)).inverse())
      .new_denominators(z2p_op) * z2p_op;
  }

  // simplified version of cctbx.sgtbx.subgroups
  af::shared<sgtbx::space_group>
  get_subgroups(sgtbx::space_group const& lattice_group)
  {
    CCTBX_ASSERT(lattice_group.n_ltr() == 1); // must be primitive
    CCTBX_ASSERT(!lattice_group.is_centric());
    af::shared<sgtbx::space_group> result;
    for(std::size_t i_smx=0;i_smx<lattice_group.n_smx();i_smx++) {
      sgtbx::space_group group_i;
      group_i.expand_smx(lattice_group.smx()[i_smx]);
      for(std::size_t j_smx=i_smx;j_smx<lattice_group.n_smx();j_smx++) {
        sgtbx::space_group subgroup(group_i);
        subgroup.expand_smx(lattice_group.smx()[j_smx]);
        subgroup.make_tidy();
        if (std::find(result.begin(),result.end(),subgroup) == result.end()) {
          result.push_back(subgroup);
        }
      }
    }
    return result;
  }

  // see cctbx.sgtbx.bravais_types
  bool
  is_bravais_type_centric(int space_group_number)
  {
    static const int numbers[] = {
      /* P -1 */ 2,
      /* P 1 2/m 1 */ 10,
      /* C 1 2/m 1 */ 12,
      /* P m m m */ 47,
      /* C m m m */ 65,
      /* F m m m */ 69,
      /* I m m m */ 71,
      /* P 4/m m m */ 123,
      /* I 4/m m m */ 139,
      /* P 6/m m m */ 191,
      /* R -3 m :H */ 166,
      /* P m -3 m */ 221,
      /* I m -3 m */ 229,
      /* F m -3 m */ 225
    };
    return std::find(numbers, numbers+14, space_group_number) != numbers+14;
  }

  // simplified version of cctbx.crystal.find_best_cell
  sgtbx::change_of_basis_op
  get_change_of_basis_op_to_best_cell(
    crystal::symmetry const& input_symmetry,
    double angular_tolerance=3)
  {
    sgtbx::change_of_basis_op best_cb_op;
    crystal::symmetry best_symmetry = input_symmetry;
    if (input_symmetry.space_group().n_smx() == 2) {
      sgtbx::rot_mx_info two_fold_info(input_symmetry.space_group()(1).r());
      CCTBX_ASSERT(scitbx::fn::absolute(two_fold_info.type()) == 2);
      scitbx::vec3<int> const& ev = two_fold_info.ev();
      CCTBX_ASSERT(std::count(ev.begin(), ev.end(), 0) == 2);
      int unique_axis = static_cast<int>(
        std::find(ev.begin(), ev.end(), 1) - ev.begin());
      sgtbx::find_affine affine(input_symmetry.space_group());
      af::const_ref<sgtbx::rt_mx> affine_cb_mx = affine.cb_mx().const_ref();
      for(std::size_t i_cb_mx=0;i_cb_mx<affine_cb_mx.size();i_cb_mx++) {
        sgtbx::change_of_basis_op cb_op = sgtbx::change_of_basis_op(
          affine_cb_mx[i_cb_mx]).new_denominators(best_cb_op);
        crystal::symmetry alt_symmetry = input_symmetry.change_basis(cb_op);
        if (alt_symmetry.space_group() == input_symmetry.space_group()) {
          int cmp_result = best_symmetry.unit_cell().compare_monoclinic(
            alt_symmetry.unit_cell(), unique_axis, angular_tolerance);
          if (cmp_result > 0) {
            best_cb_op = cb_op;
            best_symmetry = alt_symmetry;
          }
        }
      }
    }
    else {
      sgtbx::space_group affine_group("P 4 3*");
      for(std::size_t i_smx=1;i_smx<affine_group.n_smx();i_smx++) {
        sgtbx::change_of_basis_op cb_op = sgtbx::change_of_basis_op(
          affine_group(i_smx)).new_denominators(best_cb_op);
        crystal::symmetry alt_symmetry = input_symmetry.change_basis(cb_op);
        if (alt_symmetry.space_group() == input_symmetry.space_group()) {
          int cmp_result = best_symmetry.unit_cell().compare_orthorhombic(
            alt_symmetry.unit_cell());
          if (cmp_result > 0) {
            best_cb_op = cb_op;
            best_symmetry = alt_symmetry;
          }
        }
      }
    }
    return best_cb_op;
  }

  void
  show_space_group_type(sgtbx::space_group_type const& space_group_type)
  {
    std::printf("%s (No. %d)",
      space_group_type.lookup_symbol().c_str(),
      space_group_type.number());
  }

  void
  show_unit_cell(uctbx::unit_cell const& unit_cell)
  {
    for(std::size_t i=0;i<6;i++) {
      std::printf("%s%.6g%s",
        (i < 1 ? "(" : " "),
        unit_cell.parameters()[i],
        (i < 5 ? "," : ")"));
    }
  }

  // iotbx.command_line.lattice_symmetry
  void
  example(crystal::symmetry const& input_symmetry, double max_delta)
  {
    std::printf("\n");
    std::printf("Input\n");
    std::printf("=====\n");
    std::printf("\n");
    std::printf("Unit cell: ");
    show_unit_cell(input_symmetry.unit_cell());
    std::printf("\n");
    std::printf("Space group: ");
    show_space_group_type(
      sgtbx::space_group_type(input_symmetry.space_group()));
    std::printf("\n");
    std::printf("\n");
    std::printf("Angular tolerance: %.3f degrees\n", max_delta);
    std::printf("\n");
    std::printf("Similar symmetries\n");
    std::printf("==================\n");
    std::printf("\n");

    // Get cell reduction operator
    sgtbx::change_of_basis_op cb_op_inp_minimum = get_cb_op_inp_minimum(
      input_symmetry);

    // New symmetry object with changed basis
    crystal::symmetry minimum_symmetry = input_symmetry.change_basis(
      cb_op_inp_minimum);

    // precomputes potential axes (slow)
    // for repeated use store and reuse the group_search instance to save time
    sgtbx::lattice_symmetry::group_search lattice_symmetry_group;

    // Get highest symmetry compatible with lattice
    sgtbx::space_group lattice_group = lattice_symmetry_group(
      minimum_symmetry.unit_cell(), max_delta);
    lattice_group.make_tidy();

    // Get list of sub-spacegroups
    af::shared<sgtbx::space_group> subgrs = get_subgroups(lattice_group);

    // Order sub-groups
    std::vector<std::size_t> sort_values;
    sort_values.reserve(subgrs.size());
    for(std::size_t i_subgr=0;i_subgr<subgrs.size();i_subgr++) {
      std::size_t order_z = subgrs[i_subgr].order_z();
      int space_group_number = sgtbx::space_group_type(subgrs[i_subgr], false)
        .number();
      CCTBX_ASSERT(1 <= space_group_number && space_group_number <= 230);
      sort_values.push_back(order_z*1000+space_group_number);
    }
    af::shared<std::size_t> perm = af::sort_permutation(
      af::make_ref(sort_values), true);

    // Loop sub-groups in sorted order
    for(std::size_t i_perm=0;i_perm<perm.size();i_perm++) {
      sgtbx::space_group const& acentric_subgroup = subgrs[perm[i_perm]];
      // Add centre of inversion to acentric lattice symmetry
      sgtbx::space_group centric_group = sgtbx::space_group(acentric_subgroup);
      centric_group.expand_inv(sgtbx::tr_vec(0,0,0));
      // Make symmetry object: unit-cell + space-group
      // The unit cell is potentially modified to be exactly compatible
      // with the space group symmetry.
      crystal::symmetry subsym(
        centric_group.average_unit_cell(minimum_symmetry.unit_cell()),
        centric_group);
      // Ignore unwanted groups
      sgtbx::space_group_type subsym_type(subsym.space_group());
      if (!is_bravais_type_centric(subsym_type.number())) {
        continue;
      }
      // Convert subgroup to reference setting
      sgtbx::change_of_basis_op cb_op_minimum_ref = subsym_type.cb_op();
      crystal::symmetry ref_subsym = subsym.change_basis(cb_op_minimum_ref);
      // Choose best setting for monoclinic and orthorhombic systems
      sgtbx::change_of_basis_op cb_op_best_cell;
      if (3 <= subsym_type.number() && subsym_type.number() < 75) {
        cb_op_best_cell = get_change_of_basis_op_to_best_cell(ref_subsym);
      }
      crystal::symmetry best_subsym = ref_subsym.change_basis(cb_op_best_cell);
      // Total basis transformation
      sgtbx::change_of_basis_op cb_op_inp_best = cb_op_best_cell
                                               * cb_op_minimum_ref
                                               * cb_op_inp_minimum;
      // Use identity change-of-basis operator if possible
      if (best_subsym.unit_cell().is_similar_to(input_symmetry.unit_cell())) {
        sgtbx::change_of_basis_op cb_op_corr = cb_op_inp_best.inverse();
        if (   best_subsym.change_basis(cb_op_corr).space_group()
            == best_subsym.space_group()) {
          cb_op_inp_best = cb_op_corr * cb_op_inp_best;
        }
      }

      std::printf("Symmetry in minimum-lengths cell: ");
      show_space_group_type(subsym_type);
      std::printf("\n");
      std::printf("      Input minimum-lengths cell: ");
      show_unit_cell(minimum_symmetry.unit_cell());
      std::printf("\n");
      std::printf("           Symmetry-adapted cell: ");
      show_unit_cell(subsym.unit_cell());
      std::printf("\n");
      sgtbx::space_group_type best_subsym_type(best_subsym.space_group());
      std::printf("            Conventional setting: ");
      show_space_group_type(best_subsym_type);
      std::printf("\n");
      std::printf("                       Unit cell: ");
      show_unit_cell(best_subsym.unit_cell());
      std::printf("\n");
      std::printf("                 Change of basis: %s\n",
        cb_op_inp_best.c().as_xyz().c_str());
      std::printf("                         Inverse: %s\n",
        cb_op_inp_best.c_inv().as_xyz().c_str());
      std::printf("      Maximal angular difference: %.3f degrees\n",
        sgtbx::lattice_symmetry::find_max_delta(
          minimum_symmetry.unit_cell(), acentric_subgroup));
      std::printf("\n");
    }
  }

}} // namespace cctbx::<anonymous>

int main()
{
  cctbx::uctbx::unit_cell unit_cell(scitbx::af::double6(12,12,12.1,89,90,92));
  cctbx::sgtbx::space_group space_group("F 1");
  double max_delta = 3;
  cctbx::example(cctbx::crystal::symmetry(unit_cell, space_group), max_delta);
  return 0;
}
