import iotbx.pdb
from cctbx import crystal
from cctbx.crystal.direct_space_asu import non_crystallographic_asu_mappings
from cctbx.sgtbx.direct_space_asu import cut_plane
from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from cctbx.development import debug_utils
from scitbx.python_utils.misc import store
from scitbx import matrix
from scitbx.python_utils.math_utils import ifloor, iceil
from libtbx.test_utils import approx_equal
from libtbx.itertbx import count
import math
import sys

def volume_vertices(rational_asu):
  result = {}
  facets = rational_asu.facets
  n_facets = len(facets)
  for i0 in xrange(0,n_facets-2):
    for i1 in xrange(i0+1,n_facets-1):
      for i2 in xrange(i1+1,n_facets):
        m = matrix.rec(facets[i0].n+facets[i1].n+facets[i2].n,(3,3))
        d = m.determinant()
        if (d != 0):
          c = m.co_factor_matrix_transposed() / d
          b = matrix.col([-facets[i0].c,-facets[i1].c,-facets[i2].c])
          vertex = c * b
          if (rational_asu.is_inside(vertex, volume_only=0001)):
            result[vertex] = 0
  return result.keys()

def first_molecule_search_symmetry(space_group_info):
  # get the allowed origin shifts
  ss = space_group_info.structure_seminvariant()
  # assert that the continous shifts are along principal directions
  # multiply the discrete origin shifts into the group
  continuous_shifts = []
  expanded_group = sgtbx.space_group(space_group_info.type().hall_symbol())
  for vm in ss.vectors_and_moduli():
    print "ss vector:", vm.v
    print "ss modulus:", vm.m
    if (vm.m == 0):
      # collect continous allowed origin shifts
      continuous_shifts.append(vm.v)
    else:
      # add discrete allowed origin shifts
      expanded_group.expand_ltr(
        sgtbx.tr_vec(vm.v, vm.m).new_denominator(expanded_group.t_den()))
  # return the expanded group along with the continuous shifts
  return sgtbx.space_group_info(group=expanded_group), continuous_shifts

def continuous_shifts_are_principal(continuous_shifts):
  for pa in continuous_shifts:
    if (not pa in ((1,0,0),(0,1,0),(0,0,1))): return 00000
  return 0001

def fill_box(point_distance, asu_tight, asu_buffer):
  assert point_distance > 0
  f = 0.5*math.sqrt(2)*point_distance
  box_min = asu_buffer.box_min(cartesian=0001)
  box_max = asu_buffer.box_max(cartesian=0001)
  box_grid = [iceil(abs(b-e)/f) for b,e in zip(box_min,box_max)]
  print "box grid:", box_grid
  box_min = matrix.col(asu_tight.box_min(cartesian=0001))
  frac = asu_tight.unit_cell().fractionalize
  sites_cart = flex.vec3_double()
  for point in flex.nested_loop(box_grid):
    if ([p%2 for p in point].count(0) % 2 == 1):
      site_cart = matrix.col(point)*f + box_min
      site_frac = frac(site_cart)
      if (asu_buffer.is_inside(site_frac)):
        sites_cart.append(site_cart)
  return sites_cart

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

def run_one(space_group_number):
  space_group_info = sgtbx.space_group_info(symbol=space_group_number)
  expanded_group, continuous_shifts = first_molecule_search_symmetry(
    space_group_info)
  if (len(continuous_shifts) == 1):
    print "look:", continuous_shifts, space_group_info
  if (len(continuous_shifts) == 2):
    print "look:", continuous_shifts, space_group_info
  return
  space_group_info.show_summary()
  print "expanded_group:",
  expanded_group.show_summary()
  print "continuous_shifts:", continuous_shifts
  assert continuous_shifts_are_principal(continuous_shifts)
  expanded_group_type = expanded_group.type()
  reference_asu = reference_table.get_asu(expanded_group_type.number())
  reference_asu.show_comprehensive_summary()
  print expanded_group_type.cb_op().inverse().c()
  asu = reference_asu.change_basis(expanded_group_type.cb_op().inverse())
  asu.show_comprehensive_summary()
  volume_vertices = asu.volume_vertices()
  box_min = list(volume_vertices[0])
  for vertex in volume_vertices[1:]:
    for i,v in zip(count(),vertex):
      box_min[i] = min(box_min[i], v)
  box_min = matrix.col(box_min)
  for continuous_shift in continuous_shifts:
    c = -(matrix.col(continuous_shift).dot(box_min))
    cut = cut_plane.cut(n=continuous_shift, c=c)
    asu.facets.append(cut)
    asu.facets.append(-cut)
  print "asu with continuous_shifts:"
  asu.show_comprehensive_summary()
  point_distance = 2
  float_asu_tight = asu.add_buffer(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    thickness=0)
  float_asu_buffer = asu.add_buffer(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    thickness=0.5*point_distance)
  print "unit_cell:", float_asu_tight.unit_cell()
  print "box_min_frac:", float_asu_tight.box_min(cartesian=00000)
  print "box_max_frac:", float_asu_tight.box_max(cartesian=00000)
  print "box_min_cart:", float_asu_tight.box_min(cartesian=0001)
  print "box_max_cart:", float_asu_tight.box_max(cartesian=0001)
  sites_cart = fill_box(point_distance, float_asu_tight, float_asu_buffer)
  if (0):
    f = open("sites_cart.pdb", "w")
    for serial,site in zip(count(1), sites_cart):
      print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
    print >> f, "END"
    f.close()
  print "len(sites_cart):", len(sites_cart)
  check_sites(sites_cart, point_distance)
  print

def ccp_fill_box(point_distance, sampling_origin, asu_buffer):
  assert point_distance > 0
  f = 0.5*math.sqrt(2)*point_distance
  print "unit_cell:", asu_buffer.unit_cell()
  for facet in asu_buffer.facets():
    print facet.n, facet.c
  print "box_min:", asu_buffer.box_min(cartesian=0001)
  print "box_max:", asu_buffer.box_max(cartesian=0001)
  box_grid = [ifloor(abs(b-e)/f) for b,e in zip(
    asu_buffer.box_min(cartesian=0001),
    asu_buffer.box_max(cartesian=0001))]
  print "box_grid:", box_grid
  print "sampling_origin:", list(sampling_origin)
  sampling_origin = matrix.col(sampling_origin)
  frac = asu_buffer.unit_cell().fractionalize
  sites_frac = flex.vec3_double()
  for point in flex.nested_loop(box_grid):
    if ([p%2 for p in point].count(0) % 2 == 1):
      site_cart = matrix.col(point)*f + sampling_origin
      site_frac = frac(site_cart)
      if (asu_buffer.is_inside(site_frac)):
        sites_frac.append(site_frac)
  print len(sites_frac)
  return sites_frac

def cubic_close_packing_sampling(crystal_symmetry, point_distance=2):
  change_of_basis_to_reference_setting \
    = crystal_symmetry.change_of_basis_op_to_reference_setting()
  reference_symmetry = crystal_symmetry.change_basis(
    change_of_basis_to_reference_setting)
  expanded_group, continuous_shifts = first_molecule_search_symmetry(
    reference_symmetry.space_group_info())
  assert continuous_shifts_are_principal(continuous_shifts)
  rational_asu = reference_symmetry.space_group_info().direct_space_asu()
  rational_asu.add_planes(
    normal_directions=continuous_shifts,
    both_directions=0001)
  sampling_origin = rational_asu.add_buffer(
    unit_cell=reference_symmetry.unit_cell(),
    thickness=0).box_min()
  float_asu_buffer = rational_asu.add_buffer(
    unit_cell=reference_symmetry.unit_cell(),
    thickness=0.5*point_distance)
  reference_sites_frac = ccp_fill_box(
    point_distance=point_distance,
    sampling_origin=sampling_origin,
    asu_buffer=float_asu_buffer)
  rt = change_of_basis_to_reference_setting.c_inv().as_double_array()
  sites_frac = rt[:9] * reference_sites_frac
  sites_frac += rt[9:]
  return crystal_symmetry.unit_cell().orthogonalization_matrix() * sites_frac

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

def hcp_fill_box(unit_cell, point_distance, rational_asu):
  assert point_distance > 0
  hex_box = hexagonal_box(
    unit_cell=unit_cell,
    vertices_frac=volume_vertices(rational_asu=rational_asu),
    point_distance=point_distance)
  box_grid = [iceil(abs(b-e))+1 for b,e in zip(hex_box.min, hex_box.max)]
  print "box_grid:", box_grid
  float_asu_buffer = rational_asu.add_buffer(
    unit_cell=unit_cell,
    thickness=0.1*point_distance)
  hex_to_frac_matrix = (
      matrix.sqr(unit_cell.fractionalization_matrix())
    * matrix.sqr(hex_box.hexagonal_cell.orthogonalization_matrix()))
  sites_frac = flex.vec3_double()
  for point in flex.nested_loop((-1,-1,-1), box_grid, 00000):
    if (point[2] % 2 == 0):
      site_hex = [point[0],point[1],point[2]*.5]
    else:
      site_hex = [point[0]+1/3.,point[1]+2/3.,point[2]*.5]
    site_frac = hex_to_frac_matrix * matrix.col(site_hex)
    if (float_asu_buffer.is_inside(site_frac)):
      sites_frac.append(site_frac)
  return sites_frac

def hexagonal_close_packing_sampling(crystal_symmetry, point_distance=2):
  cb_op_work = crystal_symmetry.change_of_basis_op_to_reference_setting()
  point_group_type = crystal_symmetry.space_group().point_group_type()
  add_cb_op = {"2": "z,x,y",
               "m": "y,z,x"}.get(point_group_type, None)
  if (add_cb_op is not None):
    cb_op_work = sgtbx.change_of_basis_op(add_cb_op) * cb_op_work
  work_symmetry = crystal_symmetry.change_basis(cb_op_work)
  expanded_group, continuous_shifts = first_molecule_search_symmetry(
    work_symmetry.space_group_info())
  assert continuous_shifts_are_principal(continuous_shifts)
  rational_asu = work_symmetry.space_group_info().direct_space_asu()
  rational_asu.add_planes(
    normal_directions=continuous_shifts,
    both_directions=0001)
  work_sites_frac = hcp_fill_box(
    unit_cell=work_symmetry.unit_cell(),
    point_distance=point_distance,
    rational_asu=rational_asu)
  rt = cb_op_work.c_inv().as_double_array()
  sites_frac = rt[:9] * work_sites_frac
  sites_frac += rt[9:]
  return crystal_symmetry.unit_cell().orthogonalization_matrix() * sites_frac

def run_call_back(flags, space_group_info):
  crystal_symmetry = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    space_group_info=space_group_info)
  point_distance = 2
  sites_cart = hexagonal_close_packing_sampling(
    crystal_symmetry=crystal_symmetry,
    point_distance=point_distance)
  check_sites(sites_cart, point_distance)
  if (1):
    f = open("sites_cart.pdb", "w")
    for serial,site in zip(count(1), sites_cart):
      print >> f, iotbx.pdb.format_atom_record(serial=serial, site=site)
    print >> f, "END"
    f.close()

def run():
  if (0):
    run_one(1)
  if (0):
    for space_group_number in xrange(1,231):
      run_one(space_group_number)
  if (1):
    debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
