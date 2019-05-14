from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx

hexagonal_sampling_generator=crystal.close_packing_hexagonal_sampling_generator

class setup_hexagonal_sampling(object):

  def __init__(self, crystal_symmetry, symmetry_flags):
    self.cb_op_original_to_sampling = crystal_symmetry \
      .change_of_basis_op_to_reference_setting()
    point_group_type = crystal_symmetry.space_group().point_group_type()
    add_cb_op = {"2": "z,x,y",
                 "m": "y,z,x"}.get(point_group_type, None)
    if (add_cb_op is not None):
      self.cb_op_original_to_sampling = sgtbx.change_of_basis_op(add_cb_op) \
                                      * self.cb_op_original_to_sampling
    sampling_symmetry = crystal_symmetry.change_basis(
      self.cb_op_original_to_sampling)
    search_symmetry = sgtbx.search_symmetry(
      flags=symmetry_flags,
      space_group_type=sampling_symmetry.space_group_info().type(),
      seminvariant=sampling_symmetry.space_group_info()
        .structure_seminvariants())
    expanded_symmetry = crystal.symmetry(
      unit_cell=sampling_symmetry.unit_cell(),
      space_group=search_symmetry.projected_subgroup())
    self.rational_asu = expanded_symmetry.space_group_info().direct_space_asu()
    self.rational_asu.add_planes(
      normal_directions=search_symmetry.continuous_shifts(),
      both_directions=True)
    self.float_asu=self.rational_asu.define_metric(
      unit_cell=expanded_symmetry.unit_cell()).as_float_asu()
    self.continuous_shift_flags=search_symmetry.continuous_shift_flags()

def hexagonal_sampling(crystal_symmetry,
                       symmetry_flags,
                       point_distance,
                       buffer_thickness=None,
                       all_twelve_neighbors=None):
  if (buffer_thickness is None): buffer_thickness = -1
  if (all_twelve_neighbors is None): all_twelve_neighbors = False
  s = setup_hexagonal_sampling(
    crystal_symmetry=crystal_symmetry,
    symmetry_flags=symmetry_flags)
  return hexagonal_sampling_generator(
    cb_op_original_to_sampling=s.cb_op_original_to_sampling,
    float_asu=s.float_asu,
    continuous_shift_flags=s.continuous_shift_flags,
    point_distance=point_distance,
    buffer_thickness=buffer_thickness,
    all_twelve_neighbors=all_twelve_neighbors)
