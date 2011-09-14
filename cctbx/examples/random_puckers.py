"""
Simple experiment: using only bond and angle restraints of
RNA sugar ring, performing geometry minimization starting with
random coordinates, how many configurations are found?
"""

from __future__ import division
import cctbx.geometry_restraints.manager
import cctbx.geometry_restraints.lbfgs
from cctbx.array_family import flex
import scitbx.math.superpose
from scitbx import matrix
from libtbx.test_utils import approx_equal, is_below_limit
from libtbx.utils import null_out
import libtbx.load_env
import math
import sys

def pentagon_sites_cart(start_vector=(0,1.5,0), normal=(0,0,1)):
  result = flex.vec3_double([(0,0,0)])
  start_vector = matrix.col(start_vector)
  prev_point = matrix.col((0,0,0))
  axis = matrix.col(normal)
  for i in xrange(4):
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle=72*i, deg=True)
    point = prev_point + r * start_vector
    result.append(point)
    prev_point = point
  result -= result.mean()
  return result

atom_names = ["C1*", "C2*", "C3*", "C4*", "O4*"]
sites_cart_3p = flex.vec3_double([
  ( 8.601, 4.966, 0.033),
  ( 9.377, 6.271, 0.220),
  (10.287, 5.917, 1.394),
  (10.617, 4.459, 1.120),
  ( 9.403, 3.919, 0.539)])
sites_cart_2p = flex.vec3_double([
  ( 8.737, 4.960, 0.436),
  ( 9.554, 6.172, 0.859),
  (10.975, 5.614, 0.795),
  (10.788, 4.167, 1.250),
  ( 9.396, 3.833, 0.982)])
sites_cart_a = flex.vec3_double([
  (-0.39261304551933057, -1.2111974846191393, -0.1918152103139639),
  (0.39219750675929338, -0.69245062624473341, -1.397143005372731),
  (0.74981053833774947, 0.71778140495030018, -0.93667092439227351),
  (-0.49160057744166064, 1.1302189390923427, -0.1625948805803917),
  (-1.2450130796799763, -0.099026532456832558, -0.017936067449894347)])
sites_cart_b = flex.vec3_double([
  (0.47201858938285846, -0.9732520102518486, 0.33234891424712004),
  (0.83025465586861957, -0.74430566630400175, -1.1363551031011303),
  (-0.52892729188995369, -0.37071088386699969, -1.7209296245715191),
  (-1.155038969928043, 0.43588047818273268, -0.59496579629060031),
  (-0.28134783946407743, 0.20324959096386705, 0.53735380630632812)])

def run(args):
  assert args in [[], ["--verbose"]]
  if (len(args) != 0):
    cout = sys.stdout
  else:
    cout = null_out()
  edge_list_bonds = [(0,1),(0,4),(1,2),(2,3),(3,4)]
  bond_list = [
    (("C1*", "C2*"), 1.529),
    (("C1*", "O4*"), 1.412),
    (("C2*", "C3*"), 1.526),
    (("C3*", "C4*"), 1.520),
    (("C4*", "O4*"), 1.449)]
  angle_list = [
    (("C1*", "C2*", "C3*"), 101.3),
    (("C2*", "C3*", "C4*"), 102.3),
    (("C3*", "C4*", "O4*"), 104.2),
    (("C4*", "O4*", "C1*"), 110.0)]
  sites_cart, geo_manager = cctbx.geometry_restraints.manager \
    .construct_non_crystallographic_conserving_bonds_and_angles(
      sites_cart=sites_cart_3p,
      edge_list_bonds=edge_list_bonds,
      edge_list_angles=[])
  for bond_atom_names,distance_ideal in bond_list:
    i,j = [atom_names.index(atom_name) for atom_name in bond_atom_names]
    bond_params = geo_manager.bond_params_table[i][j]
    assert approx_equal(bond_params.distance_ideal, distance_ideal, eps=1.e-2)
    bond_params.distance_ideal = distance_ideal
    bond_params.weight = 1/0.02**2
  assert geo_manager.angle_proxies is None
  geo_manager.angle_proxies = cctbx.geometry_restraints.shared_angle_proxy()
  for angle_atom_names,angle_ideal in angle_list:
    i_seqs = [atom_names.index(atom_name) for atom_name in angle_atom_names]
    geo_manager.angle_proxies.append(cctbx.geometry_restraints.angle_proxy(
      i_seqs=i_seqs,
      angle_ideal=angle_ideal,
      weight=1/3**2))
  geo_manager.show_sorted(
    site_labels=atom_names, sites_cart=sites_cart, f=cout)
  def lbfgs(sites_cart):
    for i_lbfgs_restart in xrange(3):
      minimized = cctbx.geometry_restraints.lbfgs.lbfgs(
        sites_cart=sites_cart,
        geometry_restraints_manager=geo_manager)
      assert is_below_limit(value=minimized.final_target_value, limit=1e-10)
    return minimized
  lbfgs(sites_cart=sites_cart_3p)
  lbfgs(sites_cart=sites_cart_2p)
  conformer_counts = [0] * 4
  sites_cart = sites_cart.deep_copy()
  mt = flex.mersenne_twister(seed=0)
  for i_trial in xrange(20):
    while True:
      for i in xrange(sites_cart.size()):
        sites_cart[i] = mt.random_double_point_on_sphere()
      try:
        lbfgs(sites_cart=sites_cart)
      except RuntimeError, e:
        if (not str(e).startswith(
              "Bond distance > max_reasonable_bond_distance: ")):
          raise
      else:
        break
    rmsd_list = flex.double()
    for reference_sites in [
          sites_cart_3p,
          sites_cart_2p,
          sites_cart_a,
          sites_cart_b]:
      sup = scitbx.math.superpose.least_squares_fit(
        reference_sites=reference_sites,
        other_sites=sites_cart)
      rmsd = reference_sites.rms_difference(sup.other_sites_best_fit())
      rmsd_list.append(rmsd)
    oline = " ".join(["%.3f" % rmsd for rmsd in rmsd_list])
    print >> cout, oline
    assert is_below_limit(min(rmsd_list), 1e-3)
    conformer_counts[flex.min_index(rmsd_list)] += 1
  print "conformer_counts:", conformer_counts
  #
  if (libtbx.env.has_module("iotbx")):
    import iotbx.pdb.hierarchy
    hierarchy = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model(id="")
    chain = iotbx.pdb.hierarchy.chain(id="A")
    model.append_chain(chain)
    hierarchy.append_model(model)
    #
    sites_cart_pentagon = pentagon_sites_cart()
    for i_stack,sites_cart in enumerate([
          sites_cart_3p,
          sites_cart_2p,
          sites_cart_a,
          sites_cart_b]):
      atom_group = iotbx.pdb.hierarchy.atom_group(resname="  U", altloc="")
      sup = scitbx.math.superpose.least_squares_fit(
        reference_sites=sites_cart_pentagon,
        other_sites=sites_cart)
      sites_cart_out = sup.other_sites_best_fit()
      for site_label,site_cart in zip(atom_names, sites_cart_out):
        atom = iotbx.pdb.hierarchy.atom()
        atom.name = " %-3s" % site_label
        atom.xyz = matrix.col(site_cart) + matrix.col((0,0,i_stack*1.5))
        atom.occ = 1
        atom.b = 20
        atom.element = " " + site_label[0]
        atom_group.append_atom(atom)
      residue_group = iotbx.pdb.hierarchy.residue_group(
        resseq="%4d" % (i_stack+1), icode=" ")
      residue_group.append_atom_group(atom_group)
      chain.append_residue_group(residue_group)
    hierarchy.atoms().reset_serial()
    pdb_str = hierarchy.as_pdb_string(append_end=True)
    file_name = "puckers.pdb"
    print "Writing file:", file_name
    open(file_name, "w").write("""\
REMARK random_puckers.py
REMARK 1 = 3'
REMARK 2 = 2'
REMARK 3 = A
REMARK 4 = B
""" + pdb_str)
  #
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
