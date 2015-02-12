from __future__ import division
import iotbx
from scitbx.array_family import flex
from mmtbx.chemical_components import get_type

from mmtbx.monomer_library import linking_setup
from mmtbx.monomer_library.linking_setup import ad_hoc_single_metal_residue_element_types
from mmtbx.monomer_library.linking_setup import ad_hoc_non_linking_elements
from mmtbx.monomer_library.linking_setup import ad_hoc_first_row

get_class = iotbx.pdb.common_residue_names_get_class
from iotbx.pdb import get_one_letter_rna_dna_name

sugar_types = ["SACCHARIDE",
               "D-SACCHARIDE",
               "L-SACCHARIDE",
               ]
amino_types = ['"L-PEPTIDE LINKING"',
               '"D-PEPTIDE LINKING"',
               "L-PEPTIDE LINKING",
               "D-PEPTIDE LINKING",
               ]
n_linking_residues = [
  "ASN",
  ]
o_linking_residues = [
  "SER",
  "THR",
  ]
standard_n_links = [
  "NAG-ASN",
  ]
standard_o_links = [
  "NAG-SER",
  "NAG-THR",
  "MAN-SER",
  "MAN-THR",
  "XYS-SER",
  "XYS-THR",
  ]
#################################################
# saccharides that have non-standard atom names #
#  in the names in the standard links           #
#################################################
not_correct_sugars = [
  "FU4",
  ]

class empty:
  def __repr__(self):
    outl = ""
    for attr in self.__dict__:
      outl += "  %s : %s\n" % (attr, getattr(self, attr))
    return outl

def _write_warning_line(s):
  print " !!! %-78s !!!" % s

def get_distance2(atom1, atom2):
  d2 = (atom1.xyz[0]-atom2.xyz[0])**2
  d2 += (atom1.xyz[1]-atom2.xyz[1])**2
  d2 += (atom1.xyz[2]-atom2.xyz[2])**2
  return d2

def get_volume(centre, atom1, atom2, atom3):
  #if 1:
  #  from elbow.chemistry.xyzClass import xyzClass
  #  abc = []
  #  abc.append(xyzClass(atom1.xyz)-xyzClass(centre.xyz))
  #  abc.append(xyzClass(atom2.xyz)-xyzClass(centre.xyz))
  #  abc.append(xyzClass(atom3.xyz)-xyzClass(centre.xyz))
  #  volume = abc[0].DotProduct(abc[1].CrossProduct(abc[2]))
  abc = []
  abc.append(flex.double(atom1.xyz)-flex.double(centre.xyz))
  abc.append(flex.double(atom2.xyz)-flex.double(centre.xyz))
  abc.append(flex.double(atom3.xyz)-flex.double(centre.xyz))
  a = flex.vec3_double(abc[1])
  b = flex.vec3_double(abc[2])
  volume = abc[0].dot(flex.double(a.cross(b)[0]))
  return volume

def is_glyco_bond(atom1, atom2, verbose=False):
  if verbose:
    print '----- is_glyco_bond -----'
    print atom1.quote()
    print atom2.quote()
    print get_type(atom1.parent().resname)
    print get_type(atom2.parent().resname)
    print sugar_types
    print get_type(atom1.parent().resname).upper()
    print get_type(atom2.parent().resname).upper()
  if get_type(atom1.parent().resname) is None: return False
  if get_type(atom2.parent().resname) is None: return False
  if not get_type(atom1.parent().resname).upper() in sugar_types:
    return False
  if not get_type(atom2.parent().resname).upper() in sugar_types:
    return False
  #
  #if atom2.parent().resname in not_correct_sugars: return False
  return True

def is_glyco_amino_bond(atom1, atom2, verbose=False):
  if verbose:
    print '----- is_glyco_amino_bond -----'
    print atom1.quote()
    print atom2.quote()
    print get_type(atom1.parent().resname)
    print get_type(atom2.parent().resname)
    print sugar_types
    print get_type(atom1.parent().resname).upper()
    print get_type(atom2.parent().resname).upper()
  if get_type(atom1.parent().resname) is None:
    return False
  if get_type(atom2.parent().resname) is None: return False
  sugars = 0
  aminos = 0
  if get_type(atom1.parent().resname).upper() in sugar_types:
    sugars+=1
  elif get_type(atom1.parent().resname).upper() in amino_types:
    aminos+=1
  if get_type(atom2.parent().resname).upper() in sugar_types:
    sugars+=1
  elif get_type(atom2.parent().resname).upper() in amino_types:
    aminos+=1
  if sugars==1 and aminos==1:
    return True
  return False

def is_n_glyco_bond(atom1, atom2):
  if get_type(atom1.parent().resname) is None: return False
  if get_type(atom2.parent().resname) is None: return False
  sugars = 0
  n_links = 0
  if get_type(atom1.parent().resname).upper() in sugar_types:
    sugars+=1
  elif atom1.parent().resname in n_linking_residues:
    n_links+=1
  if get_type(atom2.parent().resname).upper() in sugar_types:
    sugars+=1
  elif atom2.parent().resname in n_linking_residues:
    n_links+=1
  if sugars==1 and n_links==1:
    return True
  return False

def is_o_glyco_bond(atom1, atom2):
  if get_type(atom1.parent().resname) is None: return False
  if get_type(atom2.parent().resname) is None: return False
  sugars = 0
  o_links = 0
  if get_type(atom1.parent().resname).upper() in sugar_types:
    sugars+=1
  elif atom1.parent().resname in o_linking_residues:
    o_links+=1
  if get_type(atom2.parent().resname).upper() in sugar_types:
    sugars+=1
  elif atom2.parent().resname in o_linking_residues:
    o_links+=1
  if sugars==1 and o_links==1:
    return True
  return False

def get_chiral_volume(c_atom, o_atom, angles, verbose=False):
  def _sort_by_name(a1, a2):
    if a1.name<a2.name: return -1
    else: return 1
  others = []
  for angle in angles:
    for atom in angle:
      if atom.element.strip() in ["H", "D", "T"]: continue
      if atom.parent().parent().resseq!=c_atom.parent().parent().resseq: continue
      if atom.name==c_atom.name: continue
      if atom.name==o_atom.name: continue
      others.append(atom)
  others.sort(_sort_by_name)
  others.insert(0, o_atom)
  others.insert(0, c_atom)
  if len(others)!=4:
    if verbose:
      print '-'*80
      for atom in others:
        print atom.format_atom_record()
      print '-'*80
    return None
  v = get_volume(*others)
  return v

def get_hand(c_atom, o_atom, angles, verbose=False):
  v = get_chiral_volume(c_atom, o_atom, angles, verbose=verbose)
  if v<0:
    return "BETA"
  else:
    return "ALPHA"

def get_classes(atom, important_only=False, verbose=False):
  def _num_atoms_residue(atom):
    return len(atom.parent().parent().atoms())
  def _filter_for_metal(atom, class_name):
    if class_name=="common_element":
      if atom.element.strip().upper() in ad_hoc_single_metal_residue_element_types:
        return "metal"
    if class_name=="other":
      if _num_atoms_residue(atom)==1:
        if atom.element.strip().upper() in ad_hoc_single_metal_residue_element_types:
          return "metal"
    return class_name
  #
  attrs = [
    "common_saccharide", # not in get_class
    "common_water",
    "common_element",
    "common_small_molecule",
    "common_amino_acid",
    "common_rna_dna",
    "ccp4_mon_lib_rna_dna",
    "other",
    "uncommon_amino_acid",
    "unknown",
    ]
#  elif get_type(atom1.parent().resname).upper() in amino_types:
  atom_group = atom.parent()
  classes = empty()
  for attr in attrs:
    setattr(classes, attr, False)
  if verbose:
    print '    atom_group1: altloc="%s" resname="%s" class="%s"' % (
      atom_group.altloc,
      atom_group.resname,
      get_class(atom_group.resname, consider_ccp4_mon_lib_rna_dna=True),
      )
  gc = get_class(atom_group.resname, consider_ccp4_mon_lib_rna_dna=True)
  for i, attr in enumerate(attrs):
    rc = None
    if i:
      rc = gc
    else:
      if get_type(atom_group.resname) is not None:
        if get_type(atom_group.resname).upper() in sugar_types:
          rc = attr
        #elif get_type(atom_group.resname).upper() in amino_types:
        #  rc = attr
    if rc==attr:
      if important_only: return _filter_for_metal(atom, rc)
      setattr(classes, attr, True)
  return classes

# def get_closest_atoms(atom_group1,
#                       atom_group2,
#                       ignore_hydrogens=True,
#                       ignore_atom_names_in_atom_group1=None,
#                       ignore_atom_names_in_atom_group2=None,
#                       ignore_atom_name_pairs=None,
#                       ):
#   assert 0
#   if ignore_atom_names_in_atom_group1 is None:
#     ignore_atom_names_in_atom_group1 = []
#   if ignore_atom_names_in_atom_group2 is None:
#     ignore_atom_names_in_atom_group2 = []
#   if ignore_atom_name_pairs is None: ignore_atom_name_pairs=[]
#   min_d2 = 1e5
#   min_atom1 = None
#   min_atom2 = None
#   for i, atom1 in enumerate(atom_group1.atoms()):
#     if atom1.name in ignore_atom_names_in_atom_group1: continue
#     if ignore_hydrogens:
#       if atom1.element.strip() in ad_hoc_non_linking_elements: continue
#     else:
#       if atom1.element.strip() in ad_hoc_non_linking_elements[2:]: continue
#     altloc1 = atom1.parent().altloc.strip()
#     for j, atom2 in enumerate(atom_group2.atoms()):
#       if ignore_hydrogens:
#         if atom2.element.strip() in ad_hoc_non_linking_elements: continue
#       else:
#         if atom2.element.strip() in ad_hoc_non_linking_elements[2:]: continue
#       if atom2.name in ignore_atom_names_in_atom_group2: continue
#       pair = [atom1.name, atom2.name]
#       pair.sort()
#       if pair in ignore_atom_name_pairs: continue
#       altloc2 = atom2.parent().altloc.strip()
#       if altloc1 and altloc2:
#         if altloc1!=altloc2: continue
#       d2 = get_distance2(atom1, atom2)
#       if d2<min_d2:
#         min_atom1 = atom1
#         min_atom2 = atom2
#         min_d2 = d2
#   return min_atom1, min_atom2

# def get_link_atoms(atom_group1,
#                    atom_group2,
#                    bond_cutoff=2.75,
#                    ignore_hydrogens=True,
#                    ):
#   assert 0
#   bond_cutoff *= bond_cutoff
#   link_atoms = []
#   for i, atom1 in enumerate(atom_group1.atoms()):
#     if ignore_hydrogens:
#       if atom1.element.strip() in ad_hoc_non_linking_elements: continue
#     for j, atom2 in enumerate(atom_group2.atoms()):
#       if ignore_hydrogens:
#         if atom2.element.strip() in ad_hoc_non_linking_elements: continue
#       #if i>=j: continue
#       d2 = get_distance2(atom1, atom2)
#       if d2<bond_cutoff:
#         link_atoms.append([atom1, atom2])
#   return link_atoms

# def get_nonbonded(pdb_inp,
#                   pdb_hierarchy,
#                   geometry_restraints_manager,
#                   ):
#   assert 0
#   site_labels = [atom.id_str()
#      for atom in pdb_hierarchy.atoms()]
#   pair_proxies = geometry_restraints_manager.pair_proxies(
#      sites_cart=pdb_inp.xray_structure_simple().sites_cart(),
#      site_labels=site_labels,
#      )
#   sites_cart = geometry_restraints_manager.sites_cart_used_for_pair_proxies()
#   #pair_proxies.nonbonded_proxies.show_sorted(
#   #  by_value="delta",
#   #  sites_cart=sites_cart,
#   #  )
#   site_labels = [atom.id_str()
#      for atom in pdb_hierarchy.atoms()]
#   sorted_nonbonded_proxies, not_shown = pair_proxies.nonbonded_proxies.get_sorted(
#     by_value="delta",
#     sites_cart=sites_cart,
#     site_labels=site_labels,
#     #f=sio,
#     #prefix="*",
#     #max_items=0,
#     )
#   if 0:
#     pair_proxies.nonbonded_proxies.show_sorted(
#       by_value="delta",
#       sites_cart=sites_cart,
#       site_labels=site_labels,
#       #f=sio,
#       #prefix="*",
#       #max_items=0,
#       )
#   return sorted_nonbonded_proxies

def is_atom_group_pair_linked(atom_group1,
                              atom_group2,
                              mon_lib_srv,
                              ):
  #
  # look in link list for atom group links
  #
  simple_key = "%s-%s" % (
    atom_group1.resname,
    atom_group2.resname,
    )
  if simple_key in mon_lib_srv.link_link_id_dict:
    return mon_lib_srv.link_link_id_dict[simple_key], False, simple_key
  simple_key = "%s-%s" % (
    atom_group2.resname,
    atom_group1.resname,
    )
  if simple_key in mon_lib_srv.link_link_id_dict:
    return mon_lib_srv.link_link_id_dict[simple_key], True, simple_key
  return None, None, None

def is_linked_basepairs(
        atom1,
        atom2,
        rna_dna_bond_cutoff=3.4,
        rna_dna_angle_cutoff=35):
  def final_link_direction_check():
    import math
    a1p = atom1.parent().get_atom('C4')
    a2p = atom1.parent().get_atom('C5')
    a3p = atom1.parent().get_atom('C6')
    v1p = flex.vec3_double([(a1p.xyz)]) - flex.vec3_double([(a2p.xyz)])
    v2p = flex.vec3_double([(a2p.xyz)]) - flex.vec3_double([(a3p.xyz)])
    vn = v1p.cross(v2p)
    vl = flex.vec3_double([(atom2.xyz)]) - flex.vec3_double([(atom1.xyz)])
    cos_phi = vn.dot(vl)/vn.norm()/vl.norm()
    #print "cos_phi:", cos_phi[0], "phi:", math.acos(cos_phi[0])*360/math.pi, abs(cos_phi[0]) < 0.55
    # we have cosine between normal to plane group and link, and want this angle
    # to be around 90 degrees
    return 90 - math.degrees(math.acos(abs(cos_phi[0]))) < rna_dna_angle_cutoff
  def get_distance_atoms(name1, name2):
    return atom1.parent().get_atom(name1).distance(atom2.parent().get_atom(name2))

  if atom1.distance(atom2) > rna_dna_bond_cutoff:
    return None
  from bondlength_defaults import basepairs_lengths
  if atom1.parent().resname.strip()[-1] > atom2.parent().resname.strip()[-1]:
    t = atom1
    atom1 = atom2
    atom2 = t
  rn1 = get_one_letter_rna_dna_name(atom1.parent().resname)
  rn2 = get_one_letter_rna_dna_name(atom2.parent().resname)
  atom1_idstr = atom1.id_str()
  atom2_idstr = atom2.id_str()
  # Don't link two consecutive residues in the same chain
  try:
    if ((atom1_idstr[14:15] == atom2_idstr[14:15]) and
        abs(int(atom1_idstr[16:19]) - int(atom2_idstr[16:19])) < 2 ):
      return None
  except ValueError:
    pass
  if rn1 == 'T': rn1 = 'U'
  if rn2 == 'T': rn2 = 'U'
  an1 = atom1.name.strip()
  an2 = atom2.name.strip()
  # first round of filtration
  possible_classes = []
  for k,v in basepairs_lengths.iteritems():
    if v[0] == (rn1, rn2):
      arr_to_check = [(x[0],x[1]) for x in v[1:]]
      if (an1, an2) in arr_to_check or (an2, an1) in arr_to_check:
        possible_classes.append(k)
  if len(possible_classes) == 0:
    return None
  elif len(possible_classes) >= 1:
    return possible_classes[0] if final_link_direction_check() else None
  return None

def is_atom_pair_linked(atom1,
                        atom2,
                        distance=None,
                        max_bonded_cutoff=3.,
                        amino_acid_bond_cutoff=1.9,
                        rna_dna_bond_cutoff=3.4,
                        rna_dna_angle_cutoff=35,
                        inter_residue_bond_cutoff=1.99,
                        second_row_buffer=.5,
                        metal_coordination_cutoff=3.,
                        saccharide_bond_cutoff=3.,
                        sulfur_bond_cutoff=2.5,
                        other_bond_cutoff=2., # this is the ligand distance
                        use_only_bond_cutoff=False,
                        verbose=False,
                        ):
  if atom1.element.strip().upper() in ad_hoc_non_linking_elements:
    return False
  if atom2.element.strip().upper() in ad_hoc_non_linking_elements:
    return False
  skip_if_both = linking_setup.skip_if_both
  skip_if_longer = linking_setup.update_skip_if_longer(amino_acid_bond_cutoff,
                                                       rna_dna_bond_cutoff,
                                                       inter_residue_bond_cutoff,
                                                       saccharide_bond_cutoff,
                                                       metal_coordination_cutoff,
                                                       sulfur_bond_cutoff,
                                                       other_bond_cutoff,
                                                       )
  class1 = get_classes(atom1, important_only=True)
  class2 = get_classes(atom2, important_only=True)
  class1 = linking_setup.adjust_class(atom1, class1)
  class2 = linking_setup.adjust_class(atom2, class2)
  if ( linking_setup.sulfur_class(atom1, class1)=="sulfur" and
       linking_setup.sulfur_class(atom2, class2)=="sulfur" ) :
    class1 = 'sulfur'
    class2 = 'sulfur'
  lookup = [class1, class2]
  lookup.sort()
  if verbose: print 'lookup1',lookup,skip_if_both #.get(lookup, None)
  if lookup in skip_if_both: return False
  lookup = tuple(lookup)
  if verbose: print 'lookup2',lookup
  limit = skip_if_longer.get(lookup, None)
  if limit is not None:
    if ( atom1.element not in ad_hoc_first_row or
         atom2.element not in ad_hoc_first_row):
      limit += second_row_buffer**2 # not completely accurate
  if verbose: print 'limit',limit
  if distance:
    d2 = distance**2
  else:
    d2 = get_distance2(atom1, atom2)
  if verbose: print "d2",d2
  if limit is not None and limit<d2:
    if verbose: print 'limit < d2',limit,d2
    return False
  if use_only_bond_cutoff:
    if verbose: print 'use_only_bond_cutoff',use_only_bond_cutoff
    return True
  #
  if "common_rna_dna" in lookup and "common_amino_acid" in lookup:
    return True
  if "common_rna_dna" in lookup and "common_small_molecule" in lookup:
    return True
  if "common_amino_acid" in lookup and "common_small_molecule" in lookup:
    return True
  if "other" in lookup and "common_small_molecule" in lookup:
    return True
  #
  # DNA base-pairs
  #
  if lookup.count("common_rna_dna")+lookup.count("ccp4_mon_lib_rna_dna") == 2:
    link_class = is_linked_basepairs(
        atom1,
        atom2,
        rna_dna_bond_cutoff=rna_dna_bond_cutoff,
        rna_dna_angle_cutoff=rna_dna_angle_cutoff)
    if link_class is not None:
      # print "DO LINKING, class = ", link_class, atom1.id_str(), atom2.id_str()
      return True
  #
  # sulfur bridge
  #
  if verbose:
    print atom1.quote(),linking_setup.sulfur_class(atom1, class1)
    print atom2.quote(),linking_setup.sulfur_class(atom2, class2)
  if ( linking_setup.sulfur_class(atom1, class1)=="sulfur" and
       linking_setup.sulfur_class(atom2, class2)=="sulfur" ) :
    return True
  #
  # saccharides
  #
  if verbose: print 'checking common_saccharide',lookup
  if "common_saccharide" in lookup:
    limit = saccharide_bond_cutoff*saccharide_bond_cutoff
    if verbose: print 'd2,limit',d2,limit
    if "metal" in lookup:
      limit = metal_coordination_cutoff*metal_coordination_cutoff
    if d2>limit:
      return False
    return True
  #
  # metals
  #
  if "common_element" in lookup:
    if(atom1.element.strip().upper() in ad_hoc_single_metal_residue_element_types or
       atom2.element.strip().upper() in ad_hoc_single_metal_residue_element_types
       ):
      return True
    else:
      print 'lookup',lookup
      print atom1.quote()
      print atom2.quote()
      assert 0
  if "metal" in lookup: return True
  #
  #if class1=="common_element" and class2=="common_element":
  #  assert 0
  #if class1=="common_element" or class2=="common_element":
  #  return True
  #if d2>amino_acid_bond_cutoff: return False
  #
  # non-standard amino acids
  #
  if class1=="common_amino_acid" and class2=="common_amino_acid":
    if verbose:
      print "AMINO ACIDS",atom1.quote(), atom2.quote()
  #  assert 0
  if "other" in lookup:
    if verbose: print 'other returns True'
    return True
  if verbose: print 'drop through '*5
  return False

##   if class1=="common_element" or class2=="common_element":
##     return True
##   if d2>amino_acid_bond_cutoff: return False
##   if class1=="common_amino_acid" and class2=="common_amino_acid":
##     pass # rint "AMINO ACIDS",atom1.quote(), atom2.quote()
##   return False

def generate_atoms_from_atom_groups(atom_group1, atom_group2):
  for atom in atom_group1.atoms(): yield atom
  for atom in atom_group2.atoms(): yield atom

def get_bonded_from_atom_groups(atom_group1,
                                atom_group2,
                                bond_cutoff=None,
                                verbose=False,
                                ):
  rc = {}
  if bond_cutoff:
    bond_cutoff *= bond_cutoff
  for i, atom1 in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                            atom_group2)
                                                            ):
    for j, atom2 in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                              atom_group2)
                                                              ):
      if i>=j: continue
      if bond_cutoff:
        d2 = get_distance2(atom1, atom2)
        if d2<=bond_cutoff:
          rc.setdefault(atom1.i_seq, [])
          rc[atom1.i_seq].append(atom2)
          rc.setdefault(atom2.i_seq, [])
          rc[atom2.i_seq].append(atom1)
  return rc

def get_bonded(hierarchy,
               atom,
               bond_cutoff=None,
               verbose=False,
               ):
  atoms=None
  if bond_cutoff:
    bond_cutoff *= bond_cutoff
    atoms = []
  target_atom_group = atom.parent()
  target_residue_group = target_atom_group.parent()
  target_chain = target_residue_group.parent()
  target_model = target_chain.parent()
  for model in hierarchy.models():
    if model.id!=target_model.id: continue
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if chain.id!=target_chain.id: continue
      if verbose: print 'chain: "%s"' % chain.id
      for residue_group in chain.residue_groups():
        if residue_group.resseq!=target_residue_group.resseq: continue
        if verbose: print '  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode)
        yield_residue_group = False
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if atom_group.resname!=target_atom_group.resname: continue
          if verbose: print '    atom_group: altloc="%s" resname="%s"' % (
            atom_group.altloc, atom_group.resname)
          if bond_cutoff:
            for a in atom_group.atoms():
              if a.name==atom.name: continue
              d2 = get_distance2(atom, a)
              if d2<=bond_cutoff:
                atoms.append(a)
            return atoms
          else:
            min_d2 = 1000
            min_atom = None
            for a in atom_group.atoms():
              if a.name==atom.name: continue
              d2 = get_distance2(atom, a)
              if d2<min_d2:
                min_d2 = d2
                min_atom = a
            if min_atom:
              return min_atom
  return None

def get_angles_from_included_bonds(hierarchy,
                                   bonds,
                                   bond_cutoff=None,
                                   verbose=False,
                                   ):
  tmp = []
  for bond in bonds:
    for i, atom in enumerate(bond):
      rc = get_bonded(hierarchy,
                      atom,
                      bond_cutoff=bond_cutoff,
                      verbose=verbose,
                      )
      if rc:
        for rca in rc:
          if i:
            other = bond[0]
          else:
            other = bond[1]
          tmp.append([other, atom, rca])
  if verbose:
    for angle in tmp:
      for atom in angle:
        print atom.name,
      print get_distance2(angle[0], angle[1]),
      print get_distance2(angle[1], angle[2])
  return tmp

# def process_atom_groups_for_linking(pdb_hierarchy,
#                                     atom1,
#                                     atom2,
#                                     classes1,
#                                     classes2,
#                                     #bond_cutoff=2.75,
#                                     amino_acid_bond_cutoff=1.9,
#                                     rna_dna_bond_cutoff=3.4,
#                                     intra_residue_bond_cutoff=1.99,
#                                     verbose=False,
#                                     ):
#   assert 0
#   #bond_cutoff *= bond_cutoff
#   intra_residue_bond_cutoff *= intra_residue_bond_cutoff
#   atom_group1 = atom1.parent()
#   atom_group2 = atom2.parent()
#   residue_group1 = atom_group1.parent()
#   residue_group2 = atom_group2.parent()
#   if(atom1.element.upper().strip() in ad_hoc_single_metal_residue_element_types or
#      atom2.element.upper().strip() in ad_hoc_single_metal_residue_element_types):
#     if verbose: print "Returning None because of metal"
#     return None # if metal
#     link_atoms = get_link_atoms(residue_group1, residue_group2)
#     if link_atoms:
#       return process_atom_groups_for_linking_multiple_links(pdb_hierarchy,
#                                                             link_atoms,
#                                                             verbose=verbose,
#                                                             )
#     else: return None
#   else:
#     atom1, atom2 = get_closest_atoms(residue_group1, residue_group2)
#     #if get_distance2(atom1, atom2)>bond_cutoff:
#     if atom1 is None or atom2 is None: return None
#     if get_distance2(atom1, atom2)>intra_residue_bond_cutoff:
#       if verbose: print "atoms too far apart %s %s %0.1f %0.1f" % (
#         atom1.quote(),
#         atom2.quote(),
#         get_distance2(atom1, atom2),
#         intra_residue_bond_cutoff,
#         )
#       return None
#     return process_atom_groups_for_linking_single_link(
#       pdb_hierarchy,
#       atom1,
#       atom2,
#       intra_residue_bond_cutoff=intra_residue_bond_cutoff,
#       verbose=verbose,
#       )

def process_atom_groups_for_linking_single_link(pdb_hierarchy,
                                                atom1,
                                                atom2,
                                                intra_residue_bond_cutoff=1.99,
                                                verbose=False,
                                                ):
  if is_glyco_bond(atom1, atom2, verbose=verbose):
    # glyco bonds need to be in certain order
    if atom1.name.find("C")>-1:
      tmp_atom = atom1
      atom1 = atom2
      atom2 = tmp_atom

  elif is_glyco_amino_bond(atom1, atom2, verbose=verbose):
    # needs to be better using get_class???
    if atom2.name.find("C")>-1:
      tmp_atom = atom1
      atom1 = atom2
      atom2 = tmp_atom

  long_tmp_key = "%s:%s-%s:%s" % (atom1.parent().resname.strip(),
                                  atom1.name.strip(),
                                  atom2.parent().resname.strip(),
                                  atom2.name.strip(),
    )
  tmp_key = "%s-%s" % (atom1.parent().resname.strip(),
                       atom2.parent().resname.strip(),
    )

  if verbose:
    print "tmp_key %s" % tmp_key
    print "long_tmp_key %s" % long_tmp_key
    print atom1.quote()
    print atom2.quote()
    print is_n_glyco_bond(atom1, atom2)
    print is_o_glyco_bond(atom1, atom2)
    print is_glyco_bond(atom1, atom2)

  if is_n_glyco_bond(atom1, atom2):
    if tmp_key in standard_n_links:
      data_links = ""
    key = tmp_key
  elif is_o_glyco_bond(atom1, atom2):
    if tmp_key in standard_o_links:
      data_links = ""
    key = tmp_key
  elif is_glyco_bond(atom1, atom2, verbose=verbose):
    data_links = ""
    c_atom = None
    o_atom = None
    if atom1.name.find("C")>-1:   c_atom = atom1
    elif atom2.name.find("C")>-1: c_atom = atom2
    if atom1.name.find("O")>-1:   o_atom = atom1
    elif atom2.name.find("O")>-1: o_atom = atom2
    if c_atom and o_atom:
      angles = get_angles_from_included_bonds(
        pdb_hierarchy,
        [[atom1, atom2]],
        bond_cutoff=1.75, #intra_residue_bond_cutoff,
        )
      if verbose:
        print 'get_hand'
        print c_atom, o_atom, angles
      volume = get_chiral_volume(c_atom, o_atom, angles, verbose=verbose)
      hand = get_hand(c_atom, o_atom, angles, verbose=verbose) #"ALPHA"
      if hand is None:
        key = long_tmp_key
      else:
        data_link_key = "%s%s-%s" % (hand,
                                     c_atom.name.strip()[-1],
                                     o_atom.name.strip()[-1],
                                     )
        #if data_link_key in [
        #  "BETAB-B",
        #  ]: assert 0
        #cif_links = cif_links.replace(tmp_key, data_link_key)
        key = data_link_key
    else:
      print " %s" % ("!"*86)
      _write_warning_line("  Possible link ignored")
      _write_warning_line(atom1.format_atom_record())
      _write_warning_line(atom2.format_atom_record())
      _write_warning_line("  N-linked glycan : %s" % (is_n_glyco_bond(atom1, atom2)))
      _write_warning_line("  O-linked glycan : %s" % (is_o_glyco_bond(atom1, atom2)))
      _write_warning_line("  Glycan-glycan   : %s" % (is_glyco_bond(atom1, atom2)))
      if c_atom is None: _write_warning_line("  No carbon atom found")
      if o_atom is None: _write_warning_line("  No oxygen atom found")
      print " %s" % ("!"*86)
      #print " Distance", get_distance2(atom1, atom2)
      #assert 0
      #raise Sorry("Check input geometry")
      return None
  else:
    key = long_tmp_key

  pdbres_pair = []
  for atom in [atom1, atom2]:
    pdbres_pair.append(atom.id_str(pdbres=True))
  if verbose:
    print "key %s" % key
    print pdbres_pair
    print atom1.quote()
    print atom2.quote()
  return [pdbres_pair], [key], [(atom1, atom2)]

# def process_atom_groups_for_linking_multiple_links(pdb_hierarchy,
#                                                    link_atoms,
#                                                    verbose=False,
#                                                    ):
#   assert 0
#   def _quote(atom):
#     key = ""
#     for attr in ["name", "resname", "resseq", "altloc"]:
#       if getattr(atom, attr, None) is not None:
#         key += "%s_" % getattr(atom, attr).strip()
#       elif getattr(atom.parent(), attr, None) is not None:
#         key += "%s_" % getattr(atom.parent(), attr).strip()
#       elif getattr(atom.parent().parent(), attr, None) is not None:
#         key += "%s_" % getattr(atom.parent().parent(), attr).strip()
#       else:
#         assert 0
#     return key[:-1]

#   pdbres_pairs = []
#   keys = []
#   atoms = []
#   for atom1, atom2 in link_atoms:
#     key = "%s-%s" % (_quote(atom1), _quote(atom2))
#     pdbres_pair = []
#     for atom in [atom1, atom2]:
#       pdbres_pair.append(atom.id_str(pdbres=True))
#     if verbose:
#       print atom1.quote()
#       print atom2.quote()
#       print key
#     pdbres_pairs.append(pdbres_pair)
#     keys.append(key)
#     atoms.append((atom1, atom2))
#   return pdbres_pairs, keys, atoms

def print_apply(apply):
  # from libtbx.introspection import show_stack
  outl = ''
  outl += "%s" % apply.data_link
  try:
    outl += " %s" % apply.pdbres_pair
    outl += " %s" % apply.atom1.quote()
    outl += " %s" % apply.atom2.quote()
    outl += " %s" % apply.automatic
    outl += " %s" % apply.was_used
  except Exception: pass
  #show_stack()
  return outl

# class apply_cif_list(list):
#   def __repr__(self):
#     assert 0
#     outl = "CIFs"
#     for ga in self:
#       outl += "\n%s" % print_apply(ga)
#     outl += "\n"
#     outl += '_'*80
#     return outl

#   def __append__(self, item):
#     assert 0
#     print 'APPEND'*10
#     print item
#     list.__append__(self, item)

def check_valence(hierarchy, atom):
  if atom.element.strip() not in linking_setup.simple_valence: return True
  bonds = get_bonded(hierarchy, atom, bond_cutoff=1.8)
  if len(bonds)==linking_setup.simple_valence[atom.element.strip()]:
    return False
  return True
