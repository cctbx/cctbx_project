import iotbx.pdb
from cctbx import crystal
from cctbx.crystal.direct_space_asu import non_crystallographic_asu_mappings
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from cctbx.development import debug_utils
from scitbx import matrix
from scitbx.python_utils.math_utils import iceil
from libtbx.test_utils import approx_equal
from libtbx.itertbx import count
import math
import sys

class hexagonal_box:

  def __init__(self, unit_cell, vertices_frac, point_distance):
    self.hexagonal_cell = uctbx.unit_cell((
      point_distance, point_distance, point_distance*math.sqrt(8/3.),
      90, 90, 120))
    hex_matrix = matrix.sqr(self.hexagonal_cell.fractionalization_matrix()) \
               * matrix.sqr(unit_cell.orthogonalization_matrix())
    if (len(vertices_frac) == 0):
      self.min = None
      self.max = None
    else:
      vertex_hex = hex_matrix * float(matrix.col(vertices_frac[0]))
      self.min = list(vertex_hex)
      self.max = list(vertex_hex)
      for vertex_frac in vertices_frac[1:]:
        vertex_hex = hex_matrix * float(matrix.col(vertex_frac))
        for i in xrange(3):
          self.min[i] = min(self.min[i], vertex_hex[i])
          self.max[i] = max(self.max[i], vertex_hex[i])

def hcp_fill_box(unit_cell, point_distance, rational_asu,
                 continuous_shift_flags):
  assert point_distance > 0
  hex_box = hexagonal_box(
    unit_cell=unit_cell,
    vertices_frac=rational_asu.volume_vertices(),
    point_distance=point_distance)
  box_grid = [iceil(abs(b-e))+1 for b,e in zip(hex_box.min, hex_box.max)]
  box_origin = [-1,-1,-1]
  for i,f in zip(count(), continuous_shift_flags):
    if (f):
      box_grid[i] = 0
      box_origin[i] = 0
  print "box_grid:", box_grid
  float_asu_buffer = rational_asu.add_buffer(
    unit_cell=unit_cell,
    thickness=(1+1/2.*math.sqrt(3))*point_distance)
  for facet in float_asu_buffer.facets():
    print facet.n, facet.c
  hex_to_frac_matrix = (
      matrix.sqr(unit_cell.fractionalization_matrix())
    * matrix.sqr(hex_box.hexagonal_cell.orthogonalization_matrix()))
  sites_frac = flex.vec3_double()
  for point in flex.nested_loop(begin=box_origin,
                                end=box_grid,
                                open_range=00000):
    if (point[2] % 2 == 0):
      site_hex = [point[0],point[1],point[2]*.5]
    else:
      site_hex = [point[0]+1/3.,point[1]+2/3.,point[2]*.5]
    site_frac = hex_to_frac_matrix * matrix.col(site_hex)
    if (float_asu_buffer.is_inside(site_frac)):
      sites_frac.append(site_frac)
  return sites_frac

def hexagonal_close_packing_sampling(crystal_symmetry,
                                     symmetry_flags,
                                     point_distance):
  cb_op_work = crystal_symmetry.change_of_basis_op_to_reference_setting()
  point_group_type = crystal_symmetry.space_group().point_group_type()
  add_cb_op = {"2": "z,x,y",
               "m": "y,z,x"}.get(point_group_type, None)
  if (add_cb_op is not None):
    cb_op_work = sgtbx.change_of_basis_op(add_cb_op) * cb_op_work
  work_symmetry = crystal_symmetry.change_basis(cb_op_work)
  search_symmetry = sgtbx.search_symmetry(
    flags=symmetry_flags,
    space_group_type=work_symmetry.space_group_info().type(),
    seminvariant=work_symmetry.space_group_info().structure_seminvariant())
  assert search_symmetry.continuous_shifts_are_principal()
  expanded_symmetry = crystal.symmetry(
    unit_cell=work_symmetry.unit_cell(),
    space_group=search_symmetry.group())
  rational_asu = expanded_symmetry.space_group_info().direct_space_asu()
  rational_asu.add_planes(
    normal_directions=search_symmetry.continuous_shifts(),
    both_directions=0001)
  work_sites_frac = hcp_fill_box(
    unit_cell=expanded_symmetry.unit_cell(),
    point_distance=point_distance,
    rational_asu=rational_asu,
    continuous_shift_flags=search_symmetry.continuous_shift_flags())
  rt = cb_op_work.c_inv().as_double_array()
  sites_frac = rt[:9] * work_sites_frac
  sites_frac += rt[9:]
  return crystal_symmetry.unit_cell().orthogonalization_matrix() * sites_frac

def check_sites(sites_cart, point_distance):
  asu_mappings = non_crystallographic_asu_mappings(sites_cart=sites_cart)
  distance_cutoff = point_distance * math.sqrt(2) * 0.99
  simple_pair_generator = crystal.neighbors_simple_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    full_matrix=0001)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    full_matrix=0001)
  assert simple_pair_generator.count_pairs() == pair_generator.count_pairs()
  pair_generator.restart()
  neighbors = {}
  for pair in pair_generator:
    #print pair.dist_sq**(1/2.)
    assert approx_equal(pair.dist_sq, point_distance**2)
    neighbors[pair.i_seq] = neighbors.get(pair.i_seq, 0) + 1
  n_dict = {}
  for n in neighbors.values():
    n_dict[n] = n_dict.get(n, 0) + 1
  print n_dict
  if (len(neighbors) > 0):
    assert max(neighbors.values()) <= 12

def check_with_grid_tags(inp_symmetry, symmetry_flags,
                         sites_cart, point_distance):
  cb_op_inp_ref = inp_symmetry.change_of_basis_op_to_reference_setting()
  print "cb_op_inp_ref.c():", cb_op_inp_ref.c()
  ref_symmetry = inp_symmetry.change_basis(cb_op_inp_ref)
  search_symmetry = sgtbx.search_symmetry(
    flags=symmetry_flags,
    space_group_type=ref_symmetry.space_group_info().type(),
    seminvariant=ref_symmetry.space_group_info().structure_seminvariant())
  assert search_symmetry.continuous_shifts_are_principal()
  continuous_shift_flags = search_symmetry.continuous_shift_flags()
  sites_frac_inp = inp_symmetry.unit_cell().fractionalization_matrix() \
                 * sites_cart
  rt = cb_op_inp_ref.c().as_double_array()
  sites_frac_ref = rt[:9] * sites_frac_inp
  sites_frac_ref += rt[9:]
  inp_tags = inp_symmetry.gridding(
    step=point_distance/2.,
    symmetry_flags=symmetry_flags).tags()
  for point in flex.nested_loop(inp_tags.n_real()):
    if (inp_tags.tags().tag_array()[point] < 0):
      point_frac_inp = [float(n)/d for n,d in zip(point, inp_tags.n_real())]
      point_frac_ref = cb_op_inp_ref.c() * point_frac_inp
      equiv_points = sgtbx.sym_equiv_sites(
        unit_cell=ref_symmetry.unit_cell(),
        space_group=search_symmetry.group(),
        original_site=point_frac_ref,
        minimum_distance=2.e-6,
        tolerance=1.e-6)
      min_dist = sgtbx.min_sym_equiv_distance_info(
        reference_sites=equiv_points,
        others=sites_frac_ref,
        principal_continuous_allowed_origin_shift_flags
          =continuous_shift_flags).dist()
      if (min_dist > point_distance):
        print "FAIL:", inp_symmetry.space_group_info(),point_frac_ref,min_dist
        #raise AssertionError

def run_call_back(flags, space_group_info):
  crystal_symmetry = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    space_group_info=space_group_info)
  symmetry_flags=sgtbx.search_symmetry_flags(
      use_space_group_symmetry=0001,
      use_space_group_ltr=0,
      use_seminvariant=0001,
      use_normalizer_k2l=00000,
      use_normalizer_l2n=00000)
  point_distance = 2
  sites_cart = hexagonal_close_packing_sampling(
    crystal_symmetry=crystal_symmetry,
    symmetry_flags=symmetry_flags,
    point_distance=point_distance)
  check_sites(sites_cart, point_distance)
  check_with_grid_tags(
    inp_symmetry=crystal_symmetry,
    symmetry_flags=symmetry_flags,
    sites_cart=sites_cart,
    point_distance=point_distance)
  if (0):
    f = open("sites_cart.pdb", "w")
    for serial,site in zip(count(1), sites_cart):
      print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
    print >> f, "END"
    f.close()

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
