from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.hydrogen_connectivity
import sys
import time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.model
import mmtbx.utils
from mmtbx import monomer_library
from cctbx import geometry_restraints
#import iotbx.pdb
#from scitbx.array_family import flex
#from scitbx import matrix
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx.utils import null_out
from libtbx.utils import multi_out
#from stdlib import math

legend = """\

phenix.hydrogen_connectivity:
Computes the connectivity information for hydrogen atoms.

Inputs:
  - Model file in PDB format (other file types are not supported)

Usage examples:
   phenix.hydrogen_connectivity model.pdb

Output:
  List of 1-2, 1-3 and - if necessary - 1-4 partners of all
  hydrogen atoms in the input file.

  Examples:

     HE2 :   CE2 ,  CD2  CZ  , n/a
  1-2 neighbor: CE2
  1-3 neighbors: CD2, CZ
  1-4 neighbors: not available

     HH  :   OH  ,  CZ  ,  CE2  CE1
  1-2 neighbor: OH
  1-3 neighbors: CZ
  1-4 neighbors: CE2, CE1

Important:
  Hydrogen atoms need to be present in the file.
  H atoms can be generated with phenix.reduce
"""

class atom_info(object):
  def __init__(
    self,
    iseq        = None,
    dist        = None,
    dist_ideal  = None,
    angle       = None,
    angle_ideal = None):
    adopt_init_args(self, locals())

def find_second_neighbor(ap, connectivity):
  keys = connectivity.keys()
  for i_test in ap.i_seqs:
    if(i_test in keys):
      i_h = i_test
      i_x = ((connectivity[i_h])[0][0]).iseq
      bonded = [i_h, i_x]
      i_y = [x for x in ap.i_seqs if x not in bonded][0]
      neighbor = atom_info(
        iseq = i_y,
        angle_ideal = ap.angle_ideal)
      connectivity[i_h][1].append(neighbor)

# ------------------------------------------------------------------
# This function creates dictionary, H atom is key, points to a list of lists
# {H:[[a0],[a1,a2],[b1]]}
# H        --> H atom
# [a0]     --> first neighbor (1-2)
# [a1, a2] --> second neighbors (1-3)
# [b1]     --> if only one 1-3 neighbor, list also the third neighbors
# a0, a1 are objects "atom_info"
# ------------------------------------------------------------------
def determine_H_neighbors(bond_proxies, angle_proxies, xray_structure):
  hd_selection = xray_structure.hd_selection()
  sites_cart = xray_structure.sites_cart()
  scatterers = xray_structure.scatterers()
  connectivity = {}
  # loop through bond proxies to find H atom and parent atom
  for bproxy in bond_proxies:
    i_seq, j_seq = bproxy.i_seqs
    is_i_hd = hd_selection[i_seq]
    is_j_hd = hd_selection[j_seq]
    #is there more elegant test for element?
    if(not is_i_hd and not is_j_hd): continue
    elif(is_i_hd and is_j_hd):       assert 0
    else:
      if  (is_i_hd): i_h, i_x = i_seq, j_seq
      elif(is_j_hd): i_h, i_x = j_seq, i_seq
      else:
        raise Sorry("Something went wrong in bond proxies")
    parent = atom_info(
      iseq = i_x,
      dist_ideal = bproxy.distance_ideal)
    connectivity[i_h]=[[parent]]
    connectivity[i_h].append([])
  # loop through angle proxies to find second neighbors
  for ap in angle_proxies:
    find_second_neighbor(ap, connectivity)
  # This step roughly doubles computation time
  for ih in connectivity.keys():
    reduced_neighbs = []
    for atom in connectivity[ih][1]:
      if (not hd_selection[atom.iseq]):
        reduced_neighbs.append(atom)
    if (len(reduced_neighbs) == 1):
      ix = (connectivity[ih][0][0]).iseq
      iy = (reduced_neighbs[0]).iseq
      neighbors = [ix, iy]
      connectivity[ih].append([])
      for ap in angle_proxies:
        if (ix in ap.i_seqs and iy in ap.i_seqs):
          for seq in ap.i_seqs:
            if (seq not in neighbors and not hd_selection[seq]):
              neighbor = atom_info(iseq = seq)
              connectivity[ih][2].append(neighbor)
  # TODO: eliminate double conformations in list. Might be only for HA atoms,
  # but user might come up with weird splits
  return connectivity

def run(args, out=sys.stdout):
  if (len(args) == 0):
    print legend
    return
  log = multi_out()
  log.register("stdout", out)
  log_file_name = "hydrogen_connectivity.log"
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  print >> log, "phenix.hydrogen_connectivity is running..."
  print >> log, "input parameters:\n", args

  processed_args = mmtbx.utils.process_command_line_args(
    args=args, log=null_out())
  pdb_filename = processed_args.pdb_file_names[0]
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if pdb_filename is not None:
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv    = mon_lib_srv,
      ener_lib       = ener_lib,
      file_name      = pdb_filename,
      force_symmetry = True)

  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xray_structure = processed_pdb_file.xray_structure()
  if pdb_filename is not None:
    pdb_str = pdb_hierarchy.as_pdb_string()

  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints,
    normalization = False)

  bond_proxies_simple = restraints_manager.geometry.pair_proxies(sites_cart =
    xray_structure.sites_cart()).bond_proxies.simple
  angle_proxies = restraints_manager.geometry.angle_proxies

  names = list(pdb_hierarchy.atoms().extract_name())
  sites_cart = xray_structure.sites_cart()
  scatterers = xray_structure.scatterers()

  print >>log, '\nNow determining connectivity table for H atoms...'
  connectivity = determine_H_neighbors(
    bond_proxies   = bond_proxies_simple,
    angle_proxies  = angle_proxies,
    xray_structure = xray_structure)

  if(0):
    print >>log, '\nHydrogen atom connectivity list'
    for ih in connectivity:
      if(len(connectivity[ih])==3):
        string = (" ".join([names[p.iseq] for p in connectivity[ih][2]]))
      else:
        string = 'n/a'
      print >>log, names[ih],': ', names[(connectivity[ih][0][0]).iseq], \
        ',', (" ".join([names[p.iseq] for p in connectivity[ih][1]])), \
        ',', string

if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time:", round(time.time()-t0, 2)
