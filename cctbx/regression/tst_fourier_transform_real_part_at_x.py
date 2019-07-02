from __future__ import absolute_import, division, print_function
from six.moves import zip
def set_up_random_structure(space_group_info):
  from cctbx.development import random_structure
  xray_structure = random_structure.xray_structure(
    space_group_info       = space_group_info,
    elements               =("O", "P")*4,
    volume_per_atom        = 100,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = True)
  xray_structure.scattering_type_registry(table="wk1995")
  return xray_structure

def map_value_at_sites(structure):
  densities = []
  f_calc = structure.structure_factors(d_min = 1.0).f_calc()
  phenix_fft_map = f_calc.fft_map(resolution_factor=1/6)
  phenix_fft_map.apply_volume_scaling()
  map_3d = phenix_fft_map.real_map_unpadded()
  for scatterer in structure.scatterers():
    densities.append(map_3d.eight_point_interpolation(scatterer.site))
  return densities

def map_value_at_sites_calculated(structure):
  densities = []
  f_calc = structure.structure_factors(d_min = 1.0).f_calc() \
    .generate_bijvoet_mates() \
    .expand_to_p1()
  v = f_calc.unit_cell().volume()
  for scatterer in structure.scatterers():
    rho = f_calc.indices().fourier_transform_real_part_at_x(
      fourier_coeffs=f_calc.data(), x=scatterer.site)
    densities.append(rho/v)
  return densities

def run_call_back(flags, space_group_info):
  structure = set_up_random_structure(space_group_info)
  if (flags.Verbose):
    structure.scattering_type_registry().show()
  rho_at_sites_from_phenix_fft_map = map_value_at_sites(structure)
  rho_at_sites_calculated = map_value_at_sites_calculated(structure)
  for scatterer,rf,rc in zip(
        structure.scatterers(),
        rho_at_sites_from_phenix_fft_map,
        rho_at_sites_calculated):
     if (flags.Verbose):
       print(numstr(scatterer.site), "%.3f" % rf, "%.3f" % rc)
  from scitbx.array_family import flex
  corr = flex.linear_correlation(
    flex.double(rho_at_sites_from_phenix_fft_map),
    flex.double(rho_at_sites_calculated))
  assert corr.is_well_defined
  cc = corr.coefficient()
  if (flags.Verbose):
    print("Correlation coefficient:", cc)
  from libtbx.test_utils import is_above_limit
  assert is_above_limit(value=cc, limit=0.99)
  if (flags.Verbose):
    print()

def run(args):
  from scitbx.array_family import flex
  import random
  random.seed(0)
  flex.set_random_seed(0)
  from cctbx.development import debug_utils
  debug_utils.parse_options_loop_space_groups(args, run_call_back)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
