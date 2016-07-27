from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.hydrogen_connectivity
import sys
import time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.utils
from mmtbx import monomer_library
from cctbx import geometry_restraints
from scitbx.array_family import flex
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx.utils import null_out
from libtbx.utils import multi_out
from scitbx import matrix

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
    iseq            = None,
    dist            = None,
    dist_ideal      = None,
    angle           = None,
    angle_ideal     = None,
    reduced_neighbs = None,
    count_H         = None):
    adopt_init_args(self, locals())

# ----------------------------------------------------------------------------
# This function creates a dictionary, H atom is key, points to a list of lists
# {H:[a0,[a1,a2],[b1]]}
# H        --> H atom
# a0       --> first neighbor (1-2)
# [a1, a2] --> second neighbors (1-3)
# [b1]     --> if only one 1-3 neighbor, list also the 1-4 neighbors
# a0, a1 are objects "atom_info"
# ----------------------------------------------------------------------------
def determine_H_neighbors(
  geometry_restraints, bond_proxies, angle_proxies, hd_selection, sites_cart):
  fsc2=geometry_restraints.shell_sym_tables[2].full_simple_connectivity()
  fsc1=geometry_restraints.shell_sym_tables[1].full_simple_connectivity()
  #fsc0=geometry_restraints.shell_sym_tables[0].full_simple_connectivity()
  #hd_selection = xray_structure.hd_selection()
  #sites_cart = xray_structure.sites_cart()
  # Maybe there is better way to get number of atoms?
  # (needed to select angle proxy)
  n_atoms = len(sites_cart)
  connectivity = {}
  # loop through bond proxies to find H atom and parent atom
  for bproxy in bond_proxies:
    i_seq, j_seq = bproxy.i_seqs
    is_i_hd = hd_selection[i_seq]
    is_j_hd = hd_selection[j_seq]
    if(not is_i_hd and not is_j_hd): continue
    elif(is_i_hd and is_j_hd):       assert 0
    else:
      if  (is_i_hd): ih, ix = i_seq, j_seq
      elif(is_j_hd): ih, ix = j_seq, i_seq
      else:
        raise Sorry("Something went wrong in bond proxies")
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[ix])
    dist = (r0 - rh).length()
    parent = atom_info(
      iseq       = ix,
      dist_ideal = bproxy.distance_ideal,
      dist       = dist)
    connectivity[ih]=[parent]
    connectivity[ih].append([])
    # find second neighbors
    i_parent = parent.iseq
    second_neighbors = list(fsc1[ih])
    for i_second in second_neighbors:
      iselection = flex.size_t([ih,i_parent,i_second])
      ap = angle_proxies.proxy_select(
        n_seq      = n_atoms,
        iselection = iselection)
      # Check that list ap is not empty = check if angle proxy exists
      if ap:
        angle =  ap[0].angle_ideal
        neighbor = atom_info(
          iseq = i_second,
          angle_ideal = angle)
        connectivity[ih][1].append(neighbor)
    reduced_neighbs = []
    # list of reduced neigbors is needed for finding parameterization
    for atom in connectivity[ih][1]:
      if (not hd_selection[atom.iseq]):
        reduced_neighbs.append(atom)
    (connectivity[ih][0]).reduced_neighbs = reduced_neighbs
    (connectivity[ih][0]).count_H = \
      len(connectivity[ih][1]) - len(reduced_neighbs)
    # add third neighbors, if necessary
    if (len(reduced_neighbs) == 1):
      connectivity[ih].append([])
      i_second = (reduced_neighbs[0]).iseq
      third_neighbors = list(fsc2[ih])
      for i_third in third_neighbors:
        if (not hd_selection[i_third]):
          iselection = flex.size_t([i_parent,i_second,i_third])
          ap = angle_proxies.proxy_select(
            n_seq      = n_atoms,
            iselection = iselection)
          if ap:
            neighbor = atom_info(
              iseq = i_third)
            connectivity[ih][2].append(neighbor)
    if (len(reduced_neighbs) == 2 or len(reduced_neighbs) == 3):
      ix = i_parent
      iy = (reduced_neighbs[0]).iseq
      iz = (reduced_neighbs[1]).iseq
      iselection = flex.size_t([ix,iy,iz])
      (connectivity[ih][0]).angle_ideal =  angle_proxies.proxy_select(
        n_seq      = n_atoms,
        iselection = iselection)[0].angle_ideal
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

  bond_proxies_simple, asu = restraints_manager.geometry.get_all_bond_proxies(
    sites_cart = xray_structure.sites_cart())
  angle_proxies = restraints_manager.geometry.get_all_angle_proxies()

  names = list(pdb_hierarchy.atoms().extract_name())
  sites_cart = xray_structure.sites_cart()
  #scatterers = xray_structure.scatterers()
  hd_selection = xray_structure.hd_selection()
  sites_cart = xray_structure.sites_cart()

  print >>log, '\nNow determining connectivity table for H atoms...'
  connectivity = determine_H_neighbors(
    geometry_restraints   = geometry_restraints,
    bond_proxies          = bond_proxies_simple,
    angle_proxies         = angle_proxies,
    xray_structure        = xray_structure)

  if(0):
    print >>log, '\nHydrogen atom connectivity list'
    for ih in connectivity.keys():
      if(len(connectivity[ih])==3):
        string = (" ".join([names[p.iseq] for p in connectivity[ih][2]]))
      else:
        string = 'n/a'
      print >>log, names[ih],': ', names[(connectivity[ih][0]).iseq], \
        ',', (" ".join([names[p.iseq] for p in connectivity[ih][1]])), \
        ',', string


if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time:", round(time.time()-t0, 2)
