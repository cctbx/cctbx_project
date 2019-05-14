from __future__ import absolute_import, division, print_function
import cctbx.crystal.direct_space_asu
import cctbx.sgtbx.direct_space_asu.reference_table
from cctbx.sgtbx.direct_space_asu import proto
from cctbx.array_family import flex
import sys

"""
Distance calculations using C++ classes as directly as possible.

This example is rather unusual as a Python example. There are much
higher-level Python interfaces building on the C++ classes used below
(e.g. cctbx.xray.structure.show_distances()). There is no significant
runtime penalty using the Python interfaces since all numerically
intensive calculations are implemented in C++.

The C++ implementation of the direct_space_asu class is still work in
progress (cctbx/sgtbx/direct_space_asu/proto). Currently this example
uses the Python implementation. Therefore the first part of the
find_distances() function does not directly translate into C++.
The second part of find_distances() is easily converted into C++.
"""

def find_distances(unit_cell, space_group, sites_frac, distance_cutoff):
  space_group_type = space_group.type()

  # reference_asu, metric_free_asu, and asu are Python objects
  reference_asu = cctbx.sgtbx.direct_space_asu.reference_table.get_asu(
    space_group_type.number())
  metric_free_asu = reference_asu.change_basis(
    space_group_type.cb_op() # change_of_basis_op_to_reference_setting
      .inverse())
  asu = cctbx.crystal.direct_space_asu.direct_space_asu(
      asu=metric_free_asu, unit_cell=unit_cell)

  proto_asu = proto.direct_space_asu(space_group_type) # C++ type
  # todo in C++: convert proto_asu cuts to as_float_cut_plane()

  # all objects below are wrapped C++ objects
  float_asu = cctbx.crystal.direct_space_asu_float_asu(
    unit_cell=unit_cell,
    cuts=[cut.as_float_cut_plane() for cut in metric_free_asu.cuts],
    is_inside_epsilon=1.e-6)
  asu_mappings = cctbx.crystal.direct_space_asu_asu_mappings(
    space_group=space_group,
    asu=float_asu,
    buffer_thickness=distance_cutoff)
  asu_mappings.process_sites_frac(
    original_sites=sites_frac, min_distance_sym_equiv=0.5)
  pair_asu_table = cctbx.crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  for i,pair_sym_dict in enumerate(pair_sym_table):
    print("i:", i)
    for j,sym_ops in pair_sym_dict.items():
      print("  j:", j)
      for sym_op in sym_ops:
        frac_i = sites_frac[i]
        frac_j = sites_frac[j]
        frac_ji = sym_op * frac_j
        print("    %-20s %8.3f" % (
          str(sym_op), unit_cell.distance(frac_i, frac_ji)))

def run(args):
  assert len(args) == 0
  # quartz structure
  # http://cci.lbl.gov/publications/download/iucrcompcomm_jan2003.pdf
  unit_cell = cctbx.uctbx.unit_cell((5.01,5.01,5.47,90,90,120))
  space_group = cctbx.sgtbx.space_group_info(symbol="P6222").group()
  sites_frac = flex.vec3_double([(1/2.,1/2.,1/3.), (0.197,-0.197,0.83333)])
  find_distances(
    unit_cell=unit_cell,
    space_group=space_group,
    sites_frac=sites_frac,
    distance_cutoff=5)

if (__name__ == "__main__"):
  run(sys.argv[1:])
