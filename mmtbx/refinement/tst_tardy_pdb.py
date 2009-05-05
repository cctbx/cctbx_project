import mmtbx.refinement.tardy
from mmtbx.monomer_library import pdb_interpretation
import mmtbx.monomer_library.server as mon_lib_server
import iotbx.pdb
import iotbx.phil
import cctbx.geometry_restraints
from cctbx import maptbx
from cctbx.array_family import flex
from scitbx import matrix
import libtbx.phil.command_line
from libtbx.utils import format_cpu_times
from libtbx.str_utils import format_value
from libtbx import group_args
import sys, os

class potential_object(object):

  def __init__(O,
        density_map,
        geo_manager,
        reduced_geo_manager,
        nonbonded_attenuation_factor,
        real_space_gradients_delta,
        real_space_target_weight):
    O.density_map = density_map
    O.geo_manager = geo_manager
    O.reduced_geo_manager = reduced_geo_manager
    assert isinstance(
      geo_manager.nonbonded_function,
      cctbx.geometry_restraints.prolsq_repulsion_function)
    assert nonbonded_attenuation_factor > 0
    assert nonbonded_attenuation_factor <= 1
    O.custom_nonbonded_function = geo_manager.nonbonded_function \
      .customized_copy(k_rep=nonbonded_attenuation_factor)
    O.real_space_gradients_delta = real_space_gradients_delta
    O.real_space_target_weight = real_space_target_weight
    #
    O.last_sites_moved = None
    O.f = None
    O.g = None
    O.last_grms = None

  def e_pot(O, sites_moved):
    if (O.last_sites_moved is not sites_moved):
      O.last_sites_moved = sites_moved
      sites_cart = flex.vec3_double(sites_moved)
      #
      geo_energies = O.geo_manager.energies_sites(
        sites_cart=sites_cart,
        flags=cctbx.geometry_restraints.flags.flags(nonbonded=True),
        custom_nonbonded_function=O.custom_nonbonded_function,
        compute_gradients=True)
      reduced_geo_energies = O.reduced_geo_manager.energies_sites(
        sites_cart=sites_cart,
        compute_gradients=True)
      O.f = geo_energies.target + reduced_geo_energies.target
      O.g = geo_energies.gradients + reduced_geo_energies.gradients
      O.last_grms = group_args(geo=flex.mean_sq(O.g.as_double())**0.5)
      #
      rs_f = maptbx.real_space_target_simple(
        unit_cell=O.geo_manager.crystal_symmetry.unit_cell(),
        density_map=O.density_map,
        sites_cart=sites_cart)
      rs_g = maptbx.real_space_gradients_simple(
        unit_cell=O.geo_manager.crystal_symmetry.unit_cell(),
        density_map=O.density_map,
        sites_cart=sites_cart,
        delta=O.real_space_gradients_delta)
      rs_f *= -O.real_space_target_weight
      rs_g *= -O.real_space_target_weight
      O.f += rs_f
      O.g += rs_g
      O.last_grms.real = flex.mean_sq(rs_g.as_double())**0.5
      O.last_grms.real_or_xray = "real"
      #
      O.last_grms.total = flex.mean_sq(O.g.as_double())**0.5
    return O.f

  def d_e_pot_d_sites(O, sites_moved):
    O.e_pot(sites_moved=sites_moved)
    return matrix.col_list(O.g)

master_phil_str_overrides = """\
start_temperature_kelvin = 1000
final_temperature_kelvin = 300
number_of_cooling_steps = 7
number_of_time_steps = 10
time_step_pico_seconds = 0.001
"""

def run(args, callback=None):
  master_phil = iotbx.phil.parse(
    input_string=mmtbx.refinement.tardy.master_phil_str)
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  pdb_files = []
  other_files = []
  for arg in args:
    if (os.path.isfile(arg)):
      if (iotbx.pdb.is_pdb_file(file_name=arg)):
        pdb_files.append(arg)
        arg = None
      else:
        try: phil_obj = iotbx.phil.parse(file_name=arg)
        except KeyboardInterrupt: raise
        except:
          other_files.append(arg)
        else:
          phil_objects.append(phil_obj)
        arg = None
    if (arg is not None):
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except: raise Sorry("Command-line argument not recognized: %s" % arg)
      else: phil_objects.append(command_line_params)
  params = master_phil.fetch(sources=phil_objects).extract()
  master_phil.format(params).show()
  print
  if (len(pdb_files) != 0):
    print "PDB files:"
    for file_name in pdb_files:
      print " ", file_name
    print
  if (len(other_files) != 0):
    print "Other files:"
    for file_name in other_files:
      print " ", file_name
    print
  #
  assert len(pdb_files) in [1, 2]
  #
  processed_pdb_files = pdb_interpretation.run(
    args=pdb_files[-1:]+other_files,
    strict_conflict_handling=False,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    return_all_processed_pdb_files=True)
  assert len(processed_pdb_files) == 1
  print
  #
  xs = processed_pdb_files[0].xray_structure()
  geo_manager = processed_pdb_files[0].geometry_restraints_manager()
  labels = [sc.label for sc in xs.scatterers()]
  sites = matrix.col_list(xs.sites_cart())
  masses = xs.atomic_weights()
  tardy_tree = geo_manager.construct_tardy_tree(sites=sites)
  tardy_tree.show_summary(vertex_labels=labels)
  print
  reduced_geo_manager = geo_manager.reduce_for_tardy(tardy_tree=tardy_tree)
  #
  if (len(pdb_files) == 2):
    ideal_pdb_inp = iotbx.pdb.input(file_name=pdb_files[0])
    ideal_pdb_hierarchy = ideal_pdb_inp.construct_hierarchy()
    assert ideal_pdb_hierarchy.is_similar_hierarchy(
      processed_pdb_files[0].all_chain_proxies.pdb_hierarchy)
    xs.set_sites_cart(sites_cart=ideal_pdb_hierarchy.atoms().extract_xyz())
  fft_map = xs.structure_factors(d_min=3).f_calc().fft_map()
  fft_map.apply_sigma_scaling()
  #
  potential_obj = potential_object(
    density_map=fft_map.real_map(),
    geo_manager=geo_manager,
    reduced_geo_manager=reduced_geo_manager,
    nonbonded_attenuation_factor=params.nonbonded_attenuation_factor,
    real_space_gradients_delta=0.5,
    real_space_target_weight=1)
  mmtbx.refinement.tardy.action(
    labels=labels,
    sites=sites,
    masses=masses,
    tardy_tree=tardy_tree,
    potential_obj=potential_obj,
    params=params,
    callback=callback,
    log=sys.stdout)
  print
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
