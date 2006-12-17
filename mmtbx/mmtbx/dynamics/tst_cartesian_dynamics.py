import mmtbx.restraints
from mmtbx import dynamics
from mmtbx.dynamics import cartesian_dynamics
from cctbx.array_family import flex
import time, math, os
from iotbx import pdb
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
import iotbx.pdb
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation

def exercise(test_00 = True,
             test_01 = True,
             test_02 = True,
             test_03 = True):
  try: mon_lib_srv = server.server()
  except server.MonomerLibraryServerError: return
  try: ener_lib = server.ener_lib()
  except server.MonomerLibraryServerError: return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/phe.pdb", test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping exercise(): input file not available"
    return
  processed_pdb = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=pdb_file)
  structure_initial = processed_pdb.xray_structure()
  assert structure_initial.scatterers().size() == 15
  restraints_manager = mmtbx.restraints.manager(
    geometry=processed_pdb.geometry_restraints_manager())
#
# normal run:
#
  if (test_00):
    print_flag = False
    structure_ = structure_initial.deep_copy_scatterers()
    cartesian_dynamics.cartesian_dynamics(
      structure = structure_,
      restraints_manager = restraints_manager,
      temperature = 300,
      n_steps = 200,
      time_step = 0.0005,
      verbose=-1)
    rms1 = structure_initial.rms_difference(structure_)
    rms2 = structure_.rms_difference(structure_initial)
    assert rms1 == rms2
    rms = rms1
    if(print_flag):
      print "rms between structures before and after dynamics = ", rms
    array_of_distances_between_each_atom = \
         flex.sqrt(structure_.difference_vectors_cart(structure_initial).dot())
    if(print_flag):
      print
      for d in array_of_distances_between_each_atom:
        print d
    n_rms = 4.0
    selected_by_rms = (array_of_distances_between_each_atom > n_rms * rms)
    if(n_rms > 1.0):
      assert selected_by_rms.count(True) == 0
    if(print_flag):
      print "number of outliers = ", selected_by_rms.count(True)
    selected = array_of_distances_between_each_atom.select(selected_by_rms)
    if(print_flag):
      print "list of outliers : "
      for s in selected:
        print s
#
# ran at T = 0K
#
  if (test_01):
    structure_ = structure_initial.deep_copy_scatterers()
    inst = cartesian_dynamics.cartesian_dynamics(
      structure = structure_,
      restraints_manager = restraints_manager,
      temperature = 0,
      n_steps = 200,
      time_step = 0.0005,
      verbose = -1)
    assert structure_initial.rms_difference(structure_) == structure_.rms_difference(structure_initial)
    assert approx_equal(structure_.rms_difference(structure_initial),0.0,1e-6)
#
# ran at n_step = 0
#
  if (test_02):
    structure_ = structure_initial.deep_copy_scatterers()
    cartesian_dynamics.cartesian_dynamics(
      structure = structure_,
      restraints_manager = restraints_manager,
      temperature = 300,
      n_steps = 0,
      time_step = 0.0005,
      verbose = -1)
    assert structure_initial.rms_difference(structure_) == structure_.rms_difference(structure_initial)
    assert approx_equal(structure_.rms_difference(structure_initial),0.0,1e-6)
#
# normal run with real model :
#
  if (test_03):
    pdb_file = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/2ERL_noH.pdb", test=os.path.isfile)
    if (pdb_file is None):
      print "Skipping test_03: input file not available"
      test_03 = False
  if (test_03):
    processed_pdb = pdb_interpretation.process(
                           mon_lib_srv                           = mon_lib_srv,
                           ener_lib                              = ener_lib,
                           file_name                             = pdb_file)
    structure_initial = processed_pdb.xray_structure()
    restraints_manager = mmtbx.restraints.manager(
      geometry=processed_pdb.geometry_restraints_manager())
    print_flag = False
    structure_ = structure_initial.deep_copy_scatterers()
    cartesian_dynamics.cartesian_dynamics(
                     structure                   = structure_,
                     restraints_manager = restraints_manager,
                     temperature                 = 300,
                     n_steps                     = 200,
                     time_step                   = 0.0005,
                     verbose                     = -1)
    rms1 = structure_initial.rms_difference(structure_)
    rms2 = structure_.rms_difference(structure_initial)
    assert rms1 == rms2
    rms = rms1
    if(print_flag):
      print "rms between structures before and after dynamics = ", rms
    array_of_distances_between_each_atom = \
         flex.sqrt(structure_.difference_vectors_cart(structure_initial).dot())
    if(print_flag):
      print
      for d in array_of_distances_between_each_atom:
        print d
    n_rms = 5.0
    selected_by_rms = (array_of_distances_between_each_atom > n_rms * rms)
    if(n_rms > 1.0):
      assert selected_by_rms.count(True) == 0
    if(print_flag):
      print "number of outliers = ", selected_by_rms.count(True)
    assert selected_by_rms.count(True) == 0
    selected = array_of_distances_between_each_atom.select(selected_by_rms)
    if(print_flag):
      print "list of outliers : "
      for s in selected:
        print s

def run():
  exercise()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
