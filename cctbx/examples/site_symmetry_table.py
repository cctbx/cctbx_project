from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx import uctbx

def demo():
  #
  # define unit cell parameters and a compatible space group
  #
  unit_cell = uctbx.unit_cell((10,10,15,90,90,120))
  space_group = sgtbx.space_group("-P 3 2") # Hall symbol
  #
  # define a site_symmetry_table with a few sites
  #
  site_symmetry_table = sgtbx.site_symmetry_table()
  site_symmetry_table.reserve(3) # optional: optimizes memory allocation
  for site_frac in [(0,0,0), (0.5,0.5,0.5), (0.25,0.66,0.0), (0.75,0.66,0.0)]:
    site_symmetry = sgtbx.site_symmetry(
      unit_cell=unit_cell,
      space_group=space_group,
      original_site=site_frac,
      min_distance_sym_equiv=0.5,
      assert_min_distance_sym_equiv=True)
    site_symmetry_table.process(site_symmetry_ops=site_symmetry)
  #
  # there are two sets of indices:
  #   1. "i_seq" = index into the sequence of sites as passed to
  #      site_symmetry.process().
  #   2. The indices of the tabulated special_position_ops instances.
  # site_symmetry_table.indices() establishes the relation between
  # these two sets of indices:
  #   site_symmetry_table.indices().size() = number of sites processed
  #   site_symmetry_table.indices()[i_seq] = index of special_position_ops
  #     instance in the internal site_symmetry_table.table()
  assert list(site_symmetry_table.indices()) == [1, 2, 0, 0]
  #
  # table entry 0 is always the general position
  #
  assert str(site_symmetry_table.table()[0].special_op()) == "x,y,z"
  #
  # all other table entries are special positions
  #
  assert str(site_symmetry_table.table()[1].special_op()) == "0,0,0"
  assert str(site_symmetry_table.table()[2].special_op()) == "1/2,1/2,1/2"
  #
  # To obtain the special_position_ops for a certain i_seq:
  #
  for i_seq in [0,1,2]:
    print(site_symmetry_table.get(i_seq=i_seq).special_op())
  #
  # Most of the time the (many) general positions don't need a
  # special treatment, and it is much more convenient to loop
  # only over the (few) special positions. For example, to
  # define symmetry constraints for anisotropic displacement
  # parameters:
  #
  for i_seq in site_symmetry_table.special_position_indices():
    site_constraints = site_symmetry_table.get(i_seq=i_seq).site_constraints()
    adp_constraints = site_symmetry_table.get(i_seq=i_seq).adp_constraints()
  #
  # See also:
  #   cctbx/examples/site_symmetry_constraints.py
  #   cctbx/examples/adp_symmetry_constraints.py
  #   C++ reference documentation for
  #     cctbx::sgtbx::site_symmetry_ops
  #     cctbx::sgtbx::site_symmetry
  #
  print("OK")

if (__name__ == "__main__"):
  demo()
