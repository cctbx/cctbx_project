from __future__ import division
#import sys
import math
from cctbx import geometry_restraints
from libtbx import group_args
from libtbx import adopt_init_args
from libtbx.utils import Sorry
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import dihedral_angle

# Store information about atom (angles, i_seq, ...) for h_connectivity
class atom_info(object):
  def __init__(
    self,
    iseq            = None,  # atom index number
    dist            = None,  # actual distance from coordinates
    dist_ideal      = None,  # ideal distance from restraints
    angle           = None,  # actual angle from coordinates
    angle_ideal     = None,  # ideal angle from restraints
    dihedral        = None,  # actual angle from coordinates
    dihedral_ideal  = None,  # ideal angle from restraints
    count_H         = None): # number of H atoms as second neighbors
    adopt_init_args(self, locals())

# Check if two atoms are of the same element
def is_same_element(iseq1, iseq2, atoms):
  if (atoms[iseq1].element == atoms[iseq2].element):
    result = True
  else:
    result =False
  return result

# returns dictionary: key is the index number of the atom, which points to
# the occupancy. Used to obtain conformation with highest occupancy
# {i_1:occA, i_2:occB, i_3:occC}
def make_altloc_dict(atoms, index, altloc_dict_temp, neighbors, altloc):
  altloc_dict_temp[index] = atoms[index].occ
  name = atoms[index].name
  for ag in atoms[index].parent().parent().atom_groups():
    if (ag.altloc == altloc or ag.altloc == ''):
      continue
    for ag_atom in ag.atoms():
      if (ag_atom.name == name and ag_atom.i_seq in neighbors):
        altloc_dict_temp[ag_atom.i_seq] = ag_atom.occ
  return altloc_dict_temp

# ----------------------------------------------------------------------------
# This function creates a dictionary, H atom is key, points to a list of lists
# {H:[a0,[a1,a2],[h1],[b1]]}
# H        --> H atom
# a0       --> first neighbor (1-2)
# [a1, a2] --> list second non-H neighbors (1-3)
# [h1]     --> list second H neighbors (1-3)
# [b1]     --> list 1-4 neighbors (only if necessary = only one (1-3 neighbor))
# a0, a1 are objects "atom_info"
# ----------------------------------------------------------------------------
def determine_H_neighbors(geometry_restraints, pdb_hierarchy):
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  bond_proxies_simple, asu = geometry_restraints.get_all_bond_proxies(
      sites_cart = sites_cart)
  angle_proxies = geometry_restraints.get_all_angle_proxies()
  dihedral_proxies = geometry_restraints.dihedral_proxies # this should be function in GRM, like previous
  fsc2=geometry_restraints.shell_sym_tables[2].full_simple_connectivity()
  fsc1=geometry_restraints.shell_sym_tables[1].full_simple_connectivity()
  #fsc0=geometry_restraints.shell_sym_tables[0].full_simple_h_connectivity()
  # Maybe there is better way to get number of atoms?
  n_atoms = len(sites_cart)
  h_connectivity = {}
  double_H = {}
  number_h = 0
# -------------------------------------------------------------
# loop through bond proxies to find H atom and parent atom (A0)
# -------------------------------------------------------------
  for bproxy in bond_proxies_simple:
    i_seq, j_seq = bproxy.i_seqs
    is_i_hd = atoms[i_seq].element_is_hydrogen()
    is_j_hd = atoms[j_seq].element_is_hydrogen()
    if(not is_i_hd and not is_j_hd): continue
    elif(is_i_hd and is_j_hd):       assert 0
    else:
      if  (is_i_hd): ih, i_parent = i_seq, j_seq
      elif(is_j_hd): ih, i_parent = j_seq, i_seq
      else:
        raise Sorry("Something went wrong in bond proxies")
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[i_parent])
    dist = (r0 - rh).length()
    parent = atom_info(
      iseq       = i_parent,
      dist_ideal = bproxy.distance_ideal,
      dist       = dist)
    altloc_h = atoms[ih].parent().altloc
    #print 'atom:', atoms[ih].name+' ('+str(ih)+ ') residue:', \
    #  atoms[ih].parent().parent().resseq
    i_parent_altloc = atoms[i_parent].parent().altloc
#    number_h = number_h + 1
# ---------------------------------------------------------
# if entry exists already == H has two bonds
# use the first atom encountered, ignore all others
    if (ih in h_connectivity):
      double_H[ih] = [(h_connectivity[ih][0]).iseq, i_parent]
      continue
# ---------------------------------------------------------
    number_h = number_h + 1
    # h_connectivity[ih][0] --> parent atom
    h_connectivity[ih]=[parent]
    # h_connectivity[ih][1] --> list of second non-H neighbours
    h_connectivity[ih].append([])
    # h_connectivity[ih][2] --> list of second H/D neighbours
    h_connectivity[ih].append([])
    # compute total list of second neighbors
    second_neighbors = list(fsc1[ih])
    count_H = 0
    altloc_dict = {}

    # loop to find second neighbors (ignore those where angle proxies don't
    # exist, and choose same altloc for all second neighbors)
    for i_second in second_neighbors:
      iselection = flex.size_t([ih,i_parent,i_second])
      ap = angle_proxies.proxy_select(
        n_seq      = n_atoms,
        iselection = iselection)
      if (ap):
        angle_ideal = ap[0].angle_ideal
        if (i_second in altloc_dict.keys()):
          continue
        altloc_dict_temp = {}
        i_second_altloc = atoms[i_second].parent().altloc
        # make sure that all second neighbors belong to same altloc
        if (h_connectivity[ih][1] and
          atoms[(h_connectivity[ih][1])[0].iseq].parent().altloc != ''):
          iseq_previous = (h_connectivity[ih][1])[0].iseq
          overall_altloc = atoms[iseq_previous].parent().altloc
          if (i_second_altloc != '' and altloc_h == ''
            and i_second_altloc != overall_altloc):
            continue
        # if no previous altloc, chose that which has highest occ
        elif (i_second_altloc != '' and altloc_h == ''
          and i_parent_altloc != i_second_altloc):
          altloc_dict_temp = make_altloc_dict(
            atoms            = atoms,
            index            = i_second,
            altloc_dict_temp = altloc_dict_temp,
            neighbors        = second_neighbors,
            altloc           = i_second_altloc)
        if (i_second in altloc_dict_temp.keys()):
          i_second = max(altloc_dict_temp, key=lambda k: altloc_dict_temp[k])
          iselection = flex.size_t([ih,i_parent,i_second])
          ap = angle_proxies.proxy_select(
            n_seq      = n_atoms,
            iselection = iselection)
          if ap:
            angle_ideal = ap[0].angle_ideal
        altloc_dict.update(altloc_dict_temp)
        neighbor = atom_info(
          iseq = i_second,
          angle_ideal = angle_ideal)
        is_same_hd = is_same_element(
          iseq1 = ih,
          iseq2 = i_second,
          atoms = atoms)
        if ((atoms[i_second].element_is_hydrogen() and is_same_hd) or
          (atoms[i_second].element_is_hydrogen() and i_second_altloc == '')):
          h_connectivity[ih][2].append(neighbor)
          count_H = count_H + 1
        elif (not atoms[i_second].element_is_hydrogen()):
          h_connectivity[ih][1].append(neighbor)
    (h_connectivity[ih][0]).count_H = count_H
    # ---------------------------------------------------------
    # find third neighbors, if necessary
    # ---------------------------------------------------------
    # possibly, code for second and third neighbors can be transformed
    # to a function  - TODO
    if (len(h_connectivity[ih][1]) == 1):
      h_connectivity[ih].append([])
      i_second = ((h_connectivity[ih][1])[0]).iseq
      third_neighbors = list(fsc2[ih])
      altloc_dict_third = {}
      for i_third in third_neighbors:
        if (not atoms[i_third].element_is_hydrogen()):
          iselection = flex.size_t([i_parent,i_second,i_third])
          ap = angle_proxies.proxy_select(
            n_seq      = n_atoms,
            iselection = iselection)
          if ap:
            if (i_third in altloc_dict_third.keys()):
              continue
            altloc_dict_third_temp = {}
            i_third_altloc = atoms[i_third].parent().altloc
            # make sure that all third neighbors belong to the same altloc
            if (h_connectivity[ih][3] and
              atoms[(h_connectivity[ih][3])[0].iseq].parent().altloc != ''):
              iseq_previous_third = (h_connectivity[ih][3])[0].iseq
              overall_altloc_third = atoms[iseq_previous_third].parent().altloc
              if (i_third_altloc != '' and altloc_h == ''
                and i_third_altloc != overall_altloc_third):
                continue
            elif (i_third_altloc != '' and altloc_h == ''):
              altloc_dict_third_temp = make_altloc_dict(
                atoms            = atoms,
                index            = i_third,
                altloc_dict_temp = altloc_dict_third_temp,
                neighbors        = third_neighbors,
                altloc           = i_third_altloc)
            if (i_third in altloc_dict_third_temp.keys()):
              i_third = max(altloc_dict_third_temp,
                key=lambda k: altloc_dict_third_temp[k])
            altloc_dict_third.update(altloc_dict_third_temp)
            # get dihedral angle between H, A0, A1 and third neighbor(B1, B2)
            iselection_dihe = flex.size_t([ih,i_parent,i_second,i_third])
            dp = dihedral_proxies.proxy_select(
              n_seq      = n_atoms,
              iselection = iselection_dihe)
            dihedral = dihedral_angle(
              sites=[sites_cart[i_third], sites_cart[i_second],
              sites_cart[i_parent],sites_cart[ih]])
            sites_cart_dihe = sites_cart.select(iselection_dihe).deep_copy()
            if dp:
              dihedral_id = dp[0].angle_ideal
              delta = dp.deltas(sites_cart=sites_cart_dihe)[0]
              dihedral_ideal = math.degrees(dihedral) + delta
            else:
              dihedral_ideal = None
            neighbor = atom_info(
              iseq           = i_third,
              dihedral       = dihedral,
              dihedral_ideal = dihedral_ideal)
            h_connectivity[ih][3].append(neighbor)
    # get ideal angles involving parent and other non-H second neighbors
    n_sec_neigbs = len(h_connectivity[ih][1])
    if (n_sec_neigbs == 2 or n_sec_neigbs == 3):
      reduced_neighbs = h_connectivity[ih][1]
      angles = []
      ix = i_parent
      if (n_sec_neigbs == 2):
        _list = [(0,1)]
      else:
        _list = [(0,1),(1,2),(2,0)]
      for _i,_j in _list:
        iy = (reduced_neighbs[_i]).iseq
        iz = (reduced_neighbs[_j]).iseq
        iselection = flex.size_t([ix,iy,iz])
        ap = angle_proxies.proxy_select(
          n_seq      = n_atoms,
          iselection = iselection)
        if ap:
          angles.append(ap[0].angle_ideal)
        else:
          raise Sorry("Expected ideal angles are missing. \
            Second neigbs - parent atom.")
        (h_connectivity[ih][0]).angle_ideal = angles
  #
  return group_args(
    h_connectivity = h_connectivity,
    double_H       = double_H,
    number_h       = number_h)

