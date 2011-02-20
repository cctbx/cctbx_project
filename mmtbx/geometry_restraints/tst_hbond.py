
from libtbx.test_utils import approx_equal
import cStringIO

def exercise_simple () :
  from mmtbx.geometry_restraints import hbond
  import cctbx.geometry_restraints
  from scitbx.array_family import flex
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   ALA A   1       4.288  62.025  -3.678  1.00 27.96           N
ATOM      2  CA  ALA A   1       3.119  62.736  -4.184  1.00 28.04           C
ATOM      3  C   ALA A   1       2.683  63.793  -3.199  1.00 28.16           C
ATOM      4  O   ALA A   1       2.311  64.910  -3.605  1.00 28.25           O
ATOM      5  CB  ALA A   1       1.962  61.780  -4.504  1.00 27.64           C
ATOM      6  H   ALA A   1       4.150  61.050  -3.504  1.00 27.96           H
ATOM     25  N   ALA A   5       2.586  67.940  -3.614  1.00 29.64           N
ATOM     26  CA  ALA A   5       1.415  68.721  -3.231  1.00 30.71           C
ATOM     27  C   ALA A   5       1.799  69.877  -2.302  1.00 29.77           C
ATOM     28  O   ALA A   5       1.213  70.964  -2.365  1.00 30.52           O
ATOM     29  CB  ALA A   5       0.345  67.845  -2.558  1.00 30.05           C
ATOM     30  H   ALA A   5       2.547  66.981  -3.333  1.00 29.64           H
""")
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  sites_cart = hierarchy.atoms().extract_xyz()
  cache = hierarchy.atom_selection_cache()
  i_seq_1 = cache.selection("resseq 1 and name O").iselection()
  i_seq_2 = cache.selection("resseq 5 and name H").iselection()
  assert (i_seq_1.size() == i_seq_2.size() == 1)
  build_proxies = hbond.build_simple_hbond_proxies()
  # XXX note: actual bond length is approximately 2.1021
  build_proxies.add_proxy(
    i_seqs=[i_seq_1[0], i_seq_2[0]],
    distance_ideal=1.975,
    distance_cut=2.5,
    weight=1/(0.05**2))
  cctbx_bond = cctbx.geometry_restraints.bond(
    sites=[sites_cart[i_seq_1[0]], sites_cart[i_seq_2[0]]],
    distance_ideal=1.975,
    weight=1/(0.05**2))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart)
  assert (residual == cctbx_bond.residual())
  grads = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads,
    hbond_weight=0.5)
  assert approx_equal(residual, 3.22963, eps=0.00001)
  build_proxies = hbond.build_simple_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[i_seq_1[0], i_seq_2[0]],
    distance_ideal=1.975,
    distance_cut=2.05, # unrealistic, but useful for test
    weight=1/(0.05**2))
  grads = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads)
  assert (residual == 0.0) and (grads.sum() == (0.0,0.0,0.0))
  build_proxies = hbond.build_simple_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[i_seq_1[0], i_seq_2[0]],
    distance_ideal=1.975,
    distance_cut=2.08, # unrealistic, but useful for test
    weight=1/(0.05**2))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads)
  assert approx_equal(residual, 3.6074464444, eps=0.0000001)
  #print grads.sum()
#  out = cStringIO.StringIO()
#  hbond.as_pymol_dashes(build_proxies.proxies, hierarchy, out=out)
#  print out.getvalue()
#  out = cStringIO.StringIO()
#  hbond.as_refmac_restraints(build_proxies.proxies, hierarchy, out=out,
#    sigma=0.1)
#  print out.getvalue()
  print "OK"

if (__name__ == "__main__") :
  exercise_simple()
