import mmtbx.dynamics
from mmtbx.dynamics import constants
from mmtbx.dynamics import cartesian_dynamics
import mmtbx.restraints
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from cStringIO import StringIO
import random
#from mmtbx import utils
import sys, os

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def exercise_basic(verbose):
  assert abs(constants.boltzmann_constant_akma-0.001987) < 1e-6
  assert abs(constants.akma_time_as_pico_seconds-0.04889) < 1e-5
  t = mmtbx.dynamics.kinetic_energy_as_temperature(dof=5, e=1.3)
  assert approx_equal(t, 261.678053105)
  e = mmtbx.dynamics.temperature_as_kinetic_energy(dof=5, t=t)
  assert approx_equal(e, 1.3)
  #
  masses = flex.random_double(size=10) * 4 + 1
  temps = flex.double()
  for i_pass in xrange(100):
    for i in xrange(100):
      velocities = cartesian_dynamics.random_velocities(
        masses=masses, target_temperature=300)
      kt = mmtbx.dynamics.kinetic_energy_and_temperature(velocities, masses)
      temps.append(kt.temperature)
    mmm = temps.min_max_mean()
    if (verbose): mmm.show()
    if (295 < mmm.mean < 305):
      break
  else:
    raise AssertionError("Failure reaching target_temperature.")

class get_inputs(object):

  def __init__(self, mon_lib_srv, ener_lib, verbose):
    pdb_file = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/phe.pdb", test=os.path.isfile)
    if (pdb_file is None):
      self.xray_structure = None
      self.restraints_manager = None
      return
    if (verbose): log = sys.stdout
    else:         log = StringIO()
    processed_pdb = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=pdb_file,
      log=log)
    xray_structure = processed_pdb.xray_structure()
    assert xray_structure.scatterers().size() == 15
    restraints_manager = mmtbx.restraints.manager(
      geometry=processed_pdb.geometry_restraints_manager())
    self.xray_structure = xray_structure
    self.restraints_manager = restraints_manager

def exercise_00(inputs, verbose=0):
  #
  # normal run
  #
  if (inputs.xray_structure is None):
    print "Skipping exercise_00(): input file not available"
    return
  structure_ = inputs.xray_structure.deep_copy_scatterers()
  if (verbose): log = sys.stdout
  else:         log = StringIO()
  cartesian_dynamics.cartesian_dynamics(
    structure = structure_,
    restraints_manager = inputs.restraints_manager,
    temperature = 300,
    n_steps = 200,
    time_step = 0.0005,
    log = log,
    verbose = 1)
  rms1 = inputs.xray_structure.rms_difference(structure_)
  rms2 = structure_.rms_difference(inputs.xray_structure)
  assert rms1 == rms2
  rms = rms1
  if(verbose):
    print "rms between structures before and after dynamics = ", rms
  array_of_distances_between_each_atom = \
       flex.sqrt(structure_.difference_vectors_cart(
         inputs.xray_structure).dot())
  if(verbose):
    print
    for d in array_of_distances_between_each_atom:
      print d
  n_rms = 4.0
  selected_by_rms = (array_of_distances_between_each_atom > n_rms * rms)
  if(n_rms > 1.0):
    assert selected_by_rms.count(True) == 0
  if(verbose):
    print "number of outliers = ", selected_by_rms.count(True)
  selected = array_of_distances_between_each_atom.select(selected_by_rms)
  if(verbose):
    print "list of outliers : "
    for s in selected:
      print s

def exercise_01(inputs, verbose=0):
  #
  # run at T = 0K
  #
  if (inputs.xray_structure is None):
    print "Skipping exercise_01(): input file not available"
    return
  if (verbose): log = sys.stdout
  else:         log = StringIO()
  for l,v in [(None,-1), (log, verbose)]:
    structure_ = inputs.xray_structure.deep_copy_scatterers()
    inst = cartesian_dynamics.cartesian_dynamics(
      structure = structure_,
      restraints_manager = inputs.restraints_manager,
      temperature = 0,
      n_steps = 200,
      time_step = 0.0005,
      log = l,
      verbose = v)
    assert inputs.xray_structure.rms_difference(structure_) \
        == structure_.rms_difference(inputs.xray_structure)
    assert approx_equal(
      structure_.rms_difference(inputs.xray_structure), 0.0, 1e-6)

def exercise_02(inputs, verbose=0):
  #
  # run at n_step = 0
  #
  if (inputs.xray_structure is None):
    print "Skipping exercise_02(): input file not available"
    return
  structure_ = inputs.xray_structure.deep_copy_scatterers()
  if (verbose): log = sys.stdout
  else:         log = StringIO()
  cartesian_dynamics.cartesian_dynamics(
    structure = structure_,
    restraints_manager = inputs.restraints_manager,
    temperature = 300,
    n_steps = 0,
    time_step = 0.0005,
    log = log,
    verbose = 1)
  assert inputs.xray_structure.rms_difference(structure_) \
      == structure_.rms_difference(inputs.xray_structure)
  assert approx_equal(
    structure_.rms_difference(inputs.xray_structure), 0.0, 1e-6)

def exercise_03(mon_lib_srv, ener_lib, verbose=0):
  #
  # normal run with real model
  #
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2ERL_noH.pdb", test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping exercise_03: input file not available"
    return
  if (verbose): log = sys.stdout
  else:         log = StringIO()
  processed_pdb = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv = mon_lib_srv,
    ener_lib = ener_lib,
    file_name = pdb_file,
    log = log)
  xray_structure = processed_pdb.xray_structure()
  restraints_manager = mmtbx.restraints.manager(
    geometry=processed_pdb.geometry_restraints_manager())
  structure_ = xray_structure.deep_copy_scatterers()
  cartesian_dynamics.cartesian_dynamics(
    structure = structure_,
    restraints_manager = restraints_manager,
    temperature = 300,
    n_steps = 200,
    time_step = 0.0005,
    log = log,
    verbose = 1)
  rms1 = xray_structure.rms_difference(structure_)
  rms2 = structure_.rms_difference(xray_structure)
  assert rms1 == rms2
  rms = rms1
  if(verbose):
    print "rms between structures before and after dynamics = ", rms
  array_of_distances_between_each_atom = \
       flex.sqrt(structure_.difference_vectors_cart(xray_structure).dot())
  if(verbose):
    flex.histogram(
      data=array_of_distances_between_each_atom,
      n_slots=12).show(
        format_cutoffs="%6.4f")
  n_rms = 5.2
  selected_by_rms = (array_of_distances_between_each_atom > n_rms * rms)
  outlier_sc = xray_structure.scatterers().select(selected_by_rms)
  if (outlier_sc.size() != 0):
    print "number of rms outliers:", outlier_sc.size()
    outlier_d = array_of_distances_between_each_atom.select(selected_by_rms)
    for sc,d in zip(outlier_sc, outlier_d):
      print sc.label, d
    raise RuntimeError("rms outliers.")

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_basic(verbose=verbose)
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  inputs = get_inputs(
    mon_lib_srv=mon_lib_srv, ener_lib=ener_lib, verbose=verbose)
  exercise_00(inputs=inputs, verbose=verbose)
  exercise_01(inputs=inputs, verbose=verbose)
  exercise_02(inputs=inputs, verbose=verbose)
  exercise_03(mon_lib_srv=mon_lib_srv, ener_lib=ener_lib, verbose=verbose)
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
