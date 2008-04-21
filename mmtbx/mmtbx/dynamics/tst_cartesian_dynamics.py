from mmtbx.dynamics import cartesian_dynamics
import mmtbx.restraints
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from cStringIO import StringIO
import sys, os

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib = mmtbx.monomer_library.server.ener_lib()

class get_inputs(object):

  def __init__(self, verbose):
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

def exercise_03(verbose=0):
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
    print
    for d in array_of_distances_between_each_atom:
      print d
  n_rms = 5.0
  selected_by_rms = (array_of_distances_between_each_atom > n_rms * rms)
  if(n_rms > 1.0):
    assert selected_by_rms.count(True) == 0
  if(verbose):
    print "number of outliers = ", selected_by_rms.count(True)
  assert selected_by_rms.count(True) == 0
  selected = array_of_distances_between_each_atom.select(selected_by_rms)
  if(verbose):
    print "list of outliers : "
    for s in selected:
      print s

def run():
  verbose = "--verbose" in sys.argv[1:]
  inputs = get_inputs(verbose=verbose)
  exercise_00(inputs=inputs, verbose=verbose)
  exercise_01(inputs=inputs, verbose=verbose)
  exercise_02(inputs=inputs, verbose=verbose)
  exercise_03(verbose=verbose)
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
