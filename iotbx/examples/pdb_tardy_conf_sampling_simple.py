import time
import os
op = os.path

def run(args):
  time_start = time.time()
  import iotbx.pdb
  from cctbx.eltbx.distance_based_connectivity import build_edge_list
  import scitbx.rigid_body
  import scitbx.graph.tardy_tree
  from scitbx.array_family import flex
  print "Time importing extensions: %.2f" % (time.time() - time_start)
  #
  def process(file_name):
    time_start = time.time()
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_atoms = pdb_inp.atoms()
    print "Time reading pdb file: %.2f" % (time.time() - time_start)
    print "Number of atoms:", pdb_atoms.size()
    pdb_atoms.set_chemical_element_simple_if_necessary()
    sites_cart = pdb_atoms.extract_xyz()
    #
    time_start = time.time()
    edge_list = build_edge_list(
      sites_cart=sites_cart,
      elements=pdb_atoms.extract_element())
    print "Time building bond list: %.2f" % (time.time() - time_start)
    print "Number of bonds:", len(edge_list)
    #
    time_start = time.time()
    tardy_tree = scitbx.graph.tardy_tree.construct(
      sites=sites_cart,
      edge_list=edge_list)
    print "Time building tardy tree: %.2f" % (time.time() - time_start)
    #
    time_start = time.time()
    tardy_model = scitbx.rigid_body.tardy_model(
      labels=[atom.id_str() for atom in pdb_atoms],
      sites=sites_cart,
      masses=[1]*sites_cart.size(),
      tardy_tree=tardy_tree,
      potential_obj=None)
    q = tardy_model.pack_q()
    print "Time building tardy model: %.2f" % (time.time() - time_start)
    print "Number of degrees of freedom:", q.size()
    #
    mt = flex.mersenne_twister()
    time_start = time.time()
    n_conf = 10000
    for i_conf in xrange(n_conf):
      q = mt.random_double(size=q.size())
      tardy_model.unpack_q(q_packed=q)
      conf_sites_cart = tardy_model.sites_moved()
    time_diff = time.time() - time_start
    print "time / %d conf: %.2f seconds" % (n_conf, time_diff)
    print "time / conf: %.3f milli seconds" % (time_diff / n_conf * 1000)
    if (time_diff != 0):
      print "conf / second: %.2f" % (n_conf / time_diff)
  #
  if (len(args) != 0):
    for file_name in args:
      process(file_name=file_name)
  else:
    import libtbx.load_env
    file_name = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/atp.pdb",
      test=op.isfile)
    if (file_name is None):
      from libtbx.utils import Sorry
      raise Sorry("Missing command-line argument: pdb file name")
    print "Using file:", file_name
    process(file_name=file_name)
    print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
