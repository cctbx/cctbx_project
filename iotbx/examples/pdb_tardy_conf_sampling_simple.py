"""Analyze a tardy tree from a model"""
from __future__ import absolute_import, division, print_function
import math
import time
import os
from six.moves import range
op = os.path

def build_clash_detector(n_sites, bond_list, threshold):
  import scitbx.r3_utils
  result = scitbx.r3_utils.clash_detector_simple(
    n_sites=n_sites, threshold=threshold)
  from scitbx.graph import utils
  bond_sets = utils.construct_edge_sets(
    n_vertices=n_sites, edge_list=bond_list)
  def add_exclusions(edge_sets):
    for i,edge_set in enumerate(edge_sets):
      for j in edge_set:
        if (i < j):
          result.add_exclusion(i=i, j=j)
  add_exclusions(edge_sets=bond_sets)
  angle_sets = utils.bond_bending_edge_sets(edge_sets=bond_sets)
  add_exclusions(edge_sets=angle_sets)
  for i,j in utils.potential_implied_edge_list(
               edge_sets=bond_sets, bond_bending_edge_sets=angle_sets):
    result.add_exclusion(i=i, j=j)
  return result

def run(args):
  time_start = time.time()
  import iotbx.pdb
  from cctbx.crystal.distance_based_connectivity import \
    build_simple_two_way_bond_sets
  import scitbx.rigid_body
  import scitbx.graph.tardy_tree
  from scitbx.graph.utils import extract_edge_list
  from scitbx.array_family import flex
  print("Time importing extensions: %.2f" % (time.time() - time_start))
  #
  def process(file_name, clash_threshold=2.0):
    time_start = time.time()
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_atoms = pdb_inp.atoms()
    print("Time reading pdb file: %.2f" % (time.time() - time_start))
    print("Number of atoms:", pdb_atoms.size())
    pdb_atoms.set_chemical_element_simple_if_necessary()
    sites_cart = pdb_atoms.extract_xyz()
    #
    time_start = time.time()
    bond_list = extract_edge_list(edge_sets=build_simple_two_way_bond_sets(
      sites_cart=sites_cart,
      elements=pdb_atoms.extract_element()))
    print("Time building bond list: %.2f" % (time.time() - time_start))
    print("Number of bonds:", len(bond_list))
    #
    time_start = time.time()
    tardy_tree = scitbx.graph.tardy_tree.construct(
      sites=sites_cart,
      edge_list=bond_list)
    print("Time building tardy tree: %.2f" % (time.time() - time_start))
    #
    time_start = time.time()
    tardy_model = scitbx.rigid_body.tardy_model(
      labels=[atom.id_str() for atom in pdb_atoms],
      sites=sites_cart,
      masses=[1]*sites_cart.size(),
      tardy_tree=tardy_tree,
      potential_obj=None)
    q_size_each_joint = tardy_model.q_size_each_joint()
    q_fixed = tardy_model.pack_q()[:q_size_each_joint[0]]
    assert q_size_each_joint[1:].all_eq(1) # must all be revolute joints
    q_size_moving = q_size_each_joint.size() - 1
    print("Time building tardy model: %.2f" % (time.time() - time_start))
    print("Degrees of freedom:", q_size_moving)
    #
    mt = flex.mersenne_twister()
    two_pi = 2 * math.pi
    clash_detector = build_clash_detector(
      n_sites=sites_cart.size(),
      bond_list=bond_list,
      threshold=clash_threshold)
    time_start = time.time()
    n_conf = 10000
    n_clash_conf = 0
    for i_conf in range(n_conf):
      q = q_fixed.deep_copy()
      q.extend(mt.random_double(size=q_size_moving)*two_pi)
      tardy_model.unpack_q(q_packed=q)
      conf_sites_cart = tardy_model.sites_moved()
      if (clash_detector.has_clash(sites_cart=conf_sites_cart)):
        n_clash_conf += 1
    time_diff = time.time() - time_start
    print("time / %d conf: %.2f seconds" % (n_conf, time_diff))
    print("time / conf: %.3f milli seconds" % (time_diff / n_conf * 1000))
    if (time_diff != 0):
      print("conf / second: %.2f" % (n_conf / time_diff))
    print("Fraction of conformations with clashes: %d / %d = %.2f %%" % (
      n_clash_conf, n_conf, 100. * n_clash_conf / n_conf))
  #
  if (len(args) != 0):
    for file_name in args:
      process(file_name=file_name) # PDB OK
  else:
    import libtbx.load_env
    if not libtbx.env.has_module("phenix_regression"):
      print("skipping test: phenix_regression not available.")
      return
    file_name = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/atp.pdb",
      test=op.isfile)
    if (file_name is None):
      from libtbx.utils import Sorry
      raise Sorry("Missing command-line argument: pdb file name")
    print("Using file:", file_name)
    process(file_name=file_name)
    print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
