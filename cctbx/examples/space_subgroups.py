"""
Rough script, generating space-group subgroups including super cells.
Shows the Universal Hermann-Mauguin Symbol (Zwart et al., 2008) for all
subgroups.

LIMITATION of the current implementation: settings with
fractional rotation matrix elements cannot be handled.

References:
  Grosse-Kunstleve (1999). Acta Cryst. A55, 383-395.
  Zwart et al. (2008). Acta Cryst. D64, 99-107. Section 2.1
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx
from libtbx.utils import Usage
from libtbx.str_utils import show_sorted_by_counts
from libtbx import dict_with_default_0
import sys
from six.moves import range

def loop_over_super_cells(max_index, all_subgroups, subgroup):
  assert subgroup.n_ltr() == 1
  for ia in range(1,max_index+1):
    for ib in range(1,max_index+1):
      for ic in range(1,max_index+1):
        cb_op = sgtbx.change_of_basis_op("x/%d,y/%d,z/%d" % (ia,ib,ic))
        try:
          scsubgroup = subgroup.change_basis(cb_op=cb_op)
        except RuntimeError as e:
          if (str(e).endswith(
                "Unsuitable value for rational rotation matrix.")):
            all_subgroups["incompatible_rotation_denominator"] += 1
          elif (str(e).endswith(
                "Unsuitable value for rational translation vector.")):
            all_subgroups["incompatible_translation_denominator"] += 1
          else:
            raise RuntimeError
        else:
          def remove_lattice_translations(g):
            result = sgtbx.space_group(
              hall_symbol="P1", t_den=subgroup.t_den())
            for i_inv in range(g.f_inv()):
              for i_smx in range(g.n_smx()):
                result.expand_smx(g(0, i_inv, i_smx))
            return result
          subsubgroup = remove_lattice_translations(scsubgroup)
          uhm = sgtbx.space_group_type(group=subsubgroup) \
            .universal_hermann_mauguin_symbol()
          all_subgroups[uhm] += 1

def run(args):
  if (len(args) != 2):
    raise Usage("""\
cctbx.python space_subgroups.py max_index space_group_symbol
  Example: cctbx.python space_subgroups.py 2 P41212""")
  #
  max_index = int(args[0])
  print("max_index:", max_index)
  assert max_index >= 1
  print()
  space_group_t_den = 144
  sginfo = sgtbx.space_group_info(
    symbol=args[1], space_group_t_den=space_group_t_den)
  sginfo.show_summary()
  print()
  cb_op_to_p = sginfo.change_of_basis_op_to_primitive_setting()
  sginfo_p = sginfo.change_basis(cb_op=cb_op_to_p)
  if (sginfo_p.group() != sginfo.group()):
    print("Primitive setting:")
    sginfo_p.show_summary()
    print()
  #
  all_subgroups = dict_with_default_0()
  sg_p = sginfo_p.group()
  sg_p_a = sg_p.build_derived_acentric_group()
  if (sg_p.is_centric()):
    inv_mx = sg_p(0, 1, 0).t()
  else:
    inv_mx = None
  for symx1 in sg_p_a:
    subgr1 = sgtbx.space_group(hall_symbol="P1", t_den=space_group_t_den)
    subgr1.expand_smx(symx1)
    for symx2 in sg_p_a:
      subgr2 = sgtbx.space_group(subgr1)
      subgr2.expand_smx(symx2)
      loop_over_super_cells(
        max_index=max_index, all_subgroups=all_subgroups, subgroup=subgr2)
      if (inv_mx is not None):
        subgr3 = sgtbx.space_group(subgr2)
        subgr3.expand_inv(inv_mx)
        loop_over_super_cells(
          max_index=max_index, all_subgroups=all_subgroups, subgroup=subgr3)
  #
  show_sorted_by_counts(label_count_pairs=list(all_subgroups.items()))

if (__name__ == "__main__"):
  run(sys.argv[1:])
