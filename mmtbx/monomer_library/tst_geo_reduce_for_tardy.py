from __future__ import division
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.test_utils import approx_equal
import libtbx.load_env
from cStringIO import StringIO
import sys, os

def exercise_geo_reduce_for_tardy(
      mon_lib_srv,
      ener_lib,
      file_name,
      expected_bond_counts,
      expected_dihedral_counts):
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/tardy_action/"+file_name,
    test=os.path.isfile)
  if (file_path is None):
    print 'Skipping exercise_geo_reduce_for_tardy("%s"):' \
      ' input file not available' % file_name
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_path,
    log=log)
  geo = processed_pdb_file.geometry_restraints_manager()
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart_exact()
  tardy_tree = geo.construct_tardy_tree(sites_cart=sites_cart)
  reduced_geo = geo.reduce_for_tardy(tardy_tree=tardy_tree)
  bond_counts = (
    geo.pair_proxies(sites_cart=sites_cart).bond_proxies.n_total(),
    reduced_geo.pair_proxies(sites_cart=sites_cart).bond_proxies.n_total())
  dihedral_counts = (
    geo.dihedral_proxies.size(),
    reduced_geo.dihedral_proxies.size())
  assert approx_equal(bond_counts, expected_bond_counts)
  assert approx_equal(dihedral_counts, expected_dihedral_counts)
  proxy_i_seqs_red = {}
  for proxy in reduced_geo.dihedral_proxies:
    proxy_i_seqs_red[proxy.i_seqs] = proxy
  assert len(proxy_i_seqs_red) == dihedral_counts[1]
  awl = list(processed_pdb_file.all_chain_proxies.pdb_hierarchy
    .atoms_with_labels())
  for proxy in geo.dihedral_proxies:
    if (not proxy_i_seqs_red.has_key(proxy.i_seqs)):
      sigma = 1/proxy.weight**0.5
      if (sigma > 10):
        assert awl[proxy.i_seqs[0]].resname == "PRO"

def run(args):
  assert len(args) == 0
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  for file_name,expected_bond_counts,expected_dihedral_counts in [
        ("gly_gly_box.pdb", (8,0), (2, 1)),
        ("lys_pro_trp_box.pdb", (33,0), (12, 7)),
        ("pro_lys_trp_box.pdb", (33,0), (12, 7)),
        ("1yjp_box.pdb", (59,0), (22, 15)),
        ("disulfides_box.pdb", (198,3), (56, 25))]:
    exercise_geo_reduce_for_tardy(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=file_name,
      expected_bond_counts=expected_bond_counts,
      expected_dihedral_counts=expected_dihedral_counts)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
