from __future__ import absolute_import, division, print_function
from operator import attrgetter
import iotbx
from scitbx.array_family import flex
from mmtbx.chemical_components import get_type
from libtbx.utils import Sorry

from mmtbx.monomer_library import linking_setup
from mmtbx.monomer_library.linking_setup import ad_hoc_single_metal_residue_element_types
from mmtbx.monomer_library.linking_setup import ad_hoc_non_linking_elements
from mmtbx.monomer_library.linking_setup import ad_hoc_non_linking_pairs
from mmtbx.monomer_library.linking_setup import ad_hoc_first_row

from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter

get_class = iotbx.pdb.common_residue_names_get_class

sugar_types = ["SACCHARIDE",
               "D-SACCHARIDE",
               "L-SACCHARIDE",
               'D-SACCHARIDE, ALPHA LINKING',
               'D-SACCHARIDE, BETA LINKING',
               'L-SACCHARIDE, ALPHA LINKING',
               'L-SACCHARIDE, BETA LINKING',
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

  def __lt__(self, other):
    if type(other)==type(''): return False
    for attr in sorted(self.__dict__):
      item1 = getattr(self, attr)
      item2 = getattr(other, attr)
      rc = item1==item2
      if rc:
        continue
      else:
        return rc
    return rc

def _write_warning_line(s):
  print(" !!! %-78s !!!" % s)

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
    print('----- is_glyco_bond -----')
    print(atom1.quote())
    print(atom2.quote())
    print(get_type(atom1.parent().resname))
    print(get_type(atom2.parent().resname))
    print(sugar_types)
    print(get_type(atom1.parent().resname).upper())
    print(get_type(atom2.parent().resname).upper())
    print('-------------------------')
  if get_type(atom1.parent().resname) is None: return False
  if get_type(atom2.parent().resname) is None: return False
  if not get_type(atom1.parent().resname).upper() in sugar_types:
    if verbose: print('False')
    return False
  if not get_type(atom2.parent().resname).upper() in sugar_types:
    if verbose: print('False')
    return False
  #
  #if atom2.parent().resname in not_correct_sugars: return False
  if verbose: print('True')
  return True

def is_glyco_amino_bond(atom1, atom2, verbose=False):
  if verbose:
    print('----- is_glyco_amino_bond -----')
    print(atom1.quote())
    print(atom2.quote())
    print(get_type(atom1.parent().resname))
    print(get_type(atom2.parent().resname))
    print(sugar_types)
    print(get_type(atom1.parent().resname).upper())
    print(get_type(atom2.parent().resname).upper())
    print('-------------------------------')
  if get_type(atom1.parent().resname) is None:
    if verbose: print('False')
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
    if verbose: print('True')
  if verbose: print('False')
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
  others = []
  for angle in angles:
    for atom in angle:
      if atom.element.strip() in ["H", "D", "T"]: continue
      if atom.parent().parent().resseq!=c_atom.parent().parent().resseq: continue
      if atom.name==c_atom.name: continue
      if atom.name==o_atom.name: continue
      others.append(atom)
  others.sort(key=attrgetter('name'))
  others.insert(0, o_atom)
  others.insert(0, c_atom)
  if len(others)!=4:
    if verbose:
      print('-'*80)
      for atom in others:
        print(atom.format_atom_record())
      print('-'*80)
    return None
  v = get_volume(*others)
  return v

def get_hand(c_atom, o_atom, angles, verbose=False):
  v = get_chiral_volume(c_atom, o_atom, angles, verbose=verbose)
  if v is not None and v < 0:
    return "BETA"
  else:
    return "ALPHA"

def get_classes(atom, verbose=False):
  def _num_atoms_residue(atom):
    return len(atom.parent().parent().atoms())
  #
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
  important_only_value = None
  def _filter_for_uncommon_amino_acid(atom, class_name):
    backbone_atoms = ['C', 'CA', 'N', 'O', 'OXT']
    backbone_bonds = [[0,1],]
    count = 0
    for a in atom.parent().parent().atoms():
      if a.name.strip() in backbone_atoms: count+=1
    if count>=4: class_name='uncommon_amino_acid'
    return class_name
  #
  attrs = [
    "common_saccharide",
    "common_water",
    "common_element",
    "common_small_molecule",
    "common_amino_acid",
    "common_rna_dna",
    "ccp4_mon_lib_rna_dna",
    "other",
    "uncommon_amino_acid",
    "unknown",
    'd_amino_acid',
    ]
  redirect = {"modified_amino_acid" : "other",
              "modified_rna_dna" : "other",
              }
  atom_group = atom.parent()
  classes = empty()
  for attr in attrs:
    setattr(classes, attr, False)
  # only consider ccp4 names if not single atom - CD/Cd
  consider_ccp4_mon_lib_rna_dna=True
  if len(atom_group.atoms())==1:
    consider_ccp4_mon_lib_rna_dna=False
  gc = get_class(atom_group.resname,
                 consider_ccp4_mon_lib_rna_dna=consider_ccp4_mon_lib_rna_dna)
  if verbose:
    print('    atom_group1: altloc="%s" resname="%s" class="%s"' % (
      atom_group.altloc,
      atom_group.resname,
      get_class(atom_group.resname,
                consider_ccp4_mon_lib_rna_dna=consider_ccp4_mon_lib_rna_dna),
      ))
  if atom_group.resname == 'UNK': # does this need more checks
    gc = 'common_amino_acid'
  gc = redirect.get(gc,gc)
  if verbose:
    print('final class', gc)
  for i, attr in enumerate(attrs):
    rc = None
    if i:
      rc = gc
    else:
      gotten_type = None
      if atom_group.resname in one_letter_given_three_letter:
        gotten_type = "L-PEPTIDE LINKING"
      elif atom_group.resname in ["HOH"]:
        gotten_type = "NON-POLYMER"
      #
      # special section for getting SOME of the carbohydrates using the get_class
      # or from the Chem Components which does not work in ccctbx only install
      #
      elif gc=='common_saccharide':
        rc = gc
      else:
        gotten_type = get_type(atom_group.resname)
      if gotten_type is not None:
        if gotten_type.upper() in sugar_types:
          rc = attr
      #
    if rc==attr:
      if not important_only_value:
        important_only_value = _filter_for_metal(atom, rc)
        setattr(classes, "important_only", important_only_value)
      setattr(classes, attr, True)
  if (classes.other):
    setattr(classes, _filter_for_uncommon_amino_acid(atom, rc), True)
  if verbose:
    print(classes)
  return classes

def is_atom_group_pair_linked(atom_group1,
                              atom_group2,
                              mon_lib_srv,
                              ):
  #
  # look in link list for atom group links
  #
  key_data = [
    ["%s-%s" % (atom_group1.resname, atom_group2.resname), False],
    ["%s-%s" % (atom_group2.resname, atom_group1.resname), True],
    ]
  for i, (simple_key, swap) in enumerate(key_data):
    if simple_key in mon_lib_srv.link_link_id_dict:
      return mon_lib_srv.link_link_id_dict[simple_key], swap, simple_key
  return None, None, None

def is_atom_metal_coordinated(lookup,
                              atom1,
                              atom2,
                              ):
  assert 'metal' in lookup
  if 'metal'==lookup:
    metal=atom1
    other=atom2
  else:
    metal=atom2
    other=atom1
  # print metal.quote(), other.quote()
  if other.element.strip()=='C': return False
  return True

def _get_cis_trans():
  return 'TRANS'

def allow_cis_trans_important(class_important_1, class_important_2):
  peptides = ['common_amino_acid', 'd_amino_acid', 'uncommon_amino_acid']
  if class_important_1 in peptides and class_important_2 in peptides:
    return True
  return False

def allow_cis_trans(classes1, classes2):
  return allow_cis_trans_important(classes1.important_only, classes2.important_only)

def is_atom_pair_linked_tuple(atom1,
                              atom2,
                              class_important_1,
                              class_important_2,
                              mon_lib_srv,
                              ):
  #
  # atom name specific links
  #
  if allow_cis_trans_important(class_important_1, class_important_2):
    if atom1.name==' N  ' and atom2.name==' C  ':
      return _get_cis_trans(), False, '?'
    elif atom1.name==' C  ' and atom2.name==' N  ':
      return _get_cis_trans(), True, '?'
  atom_group1 = atom1.parent()
  atom_group2 = atom2.parent()
  key_data = [
    # ['%s_%s-%s_%s' % (atom_group1.resname, atom1.name.strip(),
    #                   atom_group2.resname, atom2.name.strip()), False],
    # ['%s_%s-%s_%s' % (atom_group2.resname, atom2.name.strip(),
    #                   atom_group1.resname, atom1.name.strip()), True],
    ['%s_%s-%s_%s' % (atom_group1.resname, atom1.name.strip(),
                      'ANY', atom2.name.strip()), False],
    ['%s_%s-%s_%s' % (atom_group2.resname, atom2.name.strip(),
                      'ANY', atom1.name.strip()), True],
    ]
  for i, (simple_key, swap) in enumerate(key_data):
    if simple_key in mon_lib_srv.link_link_id_dict:
      return mon_lib_srv.link_link_id_dict[simple_key], swap, simple_key
  return None, None, None

def is_atom_pair_linked(atom1,
                        atom2,
                        class_important_1,
                        class_important_2,
                        distance=None,
                        skip_if_longer=None,
                        second_row_buffer=.5,
                        metal_coordination_cutoff=3.,
                        saccharide_bond_cutoff=3.,
                        use_only_bond_cutoff=False,
                        link_metals=True,
                        verbose=False,
                        ):
  #
  #  TODO: wrap metal_coordination_cutoff and saccharide_bond_cutoff
  #  into skip_if_longer too.
  #
  if atom1.element.strip().upper() in ad_hoc_non_linking_elements:
    return False
  if atom2.element.strip().upper() in ad_hoc_non_linking_elements:
    return False
  atom_pair = [atom1.element.strip().upper(), atom2.element.strip().upper()]
  atom_pair.sort()
  if atom_pair in ad_hoc_non_linking_pairs:
    return False
  skip_if_both = linking_setup.skip_if_both
  class1 = linking_setup.adjust_class(atom1, class_important_1)
  class2 = linking_setup.adjust_class(atom2, class_important_2)
  # python3
  # assert type(class1)==type(''), 'class1 of %s not singular : %s' % (atom1.quote(), class1)
  # assert type(class2)==type(''), 'class2 of %s not singular : %s' % (atom2.quote(), class2)
  if ( linking_setup.sulfur_class(atom1, class1)=="sulfur" and
       linking_setup.sulfur_class(atom2, class2)=="sulfur" ):
    class1 = 'sulfur'
    class2 = 'sulfur'
  lookup = [class1, class2]
  if verbose: print('lookup', lookup, atom1.quote(), atom2.quote())
  lookup.sort()
  if verbose: print('lookup1',lookup,skip_if_both) #.get(lookup, None)
  if lookup in skip_if_both: return False
  lookup = tuple(lookup)
  if verbose: print('lookup2',lookup)
  limit = skip_if_longer.get(lookup, None)
  if limit is not None:
    if ( atom1.element not in ad_hoc_first_row or
         atom2.element not in ad_hoc_first_row):
      limit += second_row_buffer**2 # not completely accurate
  if verbose: print('limit',limit)
  if distance:
    d2 = distance**2
  else:
    d2 = get_distance2(atom1, atom2)
  if verbose: print("d2",d2)
  if limit is not None and limit<d2:
    if verbose: print('limit < d2',limit,d2)
    return False
  if use_only_bond_cutoff:
    if verbose: print('use_only_bond_cutoff',use_only_bond_cutoff)
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
  # sulfur bridge
  #
  if verbose:
    print(atom1.quote(),linking_setup.sulfur_class(atom1, class1))
    print(atom2.quote(),linking_setup.sulfur_class(atom2, class2))
  if ( linking_setup.sulfur_class(atom1, class1)=="sulfur" and
       linking_setup.sulfur_class(atom2, class2)=="sulfur" ):
    return True
  #
  # saccharides
  #
  if verbose: print('checking common_saccharide',lookup)
  if "common_saccharide" in lookup:
    limit = saccharide_bond_cutoff*saccharide_bond_cutoff
    if "metal" in lookup:
      limit = metal_coordination_cutoff*metal_coordination_cutoff
    if verbose: print('d2,limit',d2,limit)
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
      outl = '''

Linking of a common element has failed.

  %-20s - %s
  %-20s - %s

Send details to help@phenix-online.org
      ''' % (lookup[0], atom1.quote(), lookup[1], atom2.quote())
      raise Sorry(outl)
  if "metal" in lookup:
    if not linking_setup.skip_if_non_linking(lookup, atom1, atom2):
      return False
    else:
      if link_metals is True:
        return is_atom_metal_coordinated(lookup, atom1, atom2)
      return link_metals
  #
  # amino acids
  #
  if class1=="common_amino_acid" and class2=="common_amino_acid":
    if verbose:
      print("AMINO ACIDS",atom1.quote(), atom2.quote())
    el1 = atom1.element.strip().upper()
    el2 = atom2.element.strip().upper()
    if ( (el1=='N' and el2=='C') or (el1=='C' and el2=='N')):
      return True

  #
  # D-peptide special case...
  #
  if class1=='d_amino_acid' or class2=='d_amino_acid':
    if class1=='d_amino_acid':
      d_amino_acid=atom1
    elif class2=='d_amino_acid':
      d_amino_acid=atom2
    if d_amino_acid.name.strip() in ['O']:
      if verbose: print('d_amino_acid do not link O')
      return False
  #
  # other
  #
  if lookup==("other", "other"):
    # non-standard residue to non-standard is too risky to do automatically ...
    if verbose: print('other, other returns True')
    return True
  elif "other" in lookup:
    if verbose: print('other returns True')
    return True
  if verbose: print('drop through '*5)
  return False

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
  for i, atom1 in enumerate(generate_atoms_from_atom_groups(atom_group1.parent(),
                                                            atom_group2.parent())
                                                            ):
    for j, atom2 in enumerate(generate_atoms_from_atom_groups(atom_group1.parent(),
                                                              atom_group2.parent())
                                                              ):
      if i>=j: continue
      if bond_cutoff:
        altloc1 = atom1.parent().altloc
        altloc2 = atom2.parent().altloc
        if altloc1==altloc2: pass
        elif altloc1=='' or altloc2=='': pass
        else: continue
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
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if chain.id!=target_chain.id: continue
      if verbose: print('chain: "%s"' % chain.id)
      for residue_group in chain.residue_groups():
        if residue_group.resseq!=target_residue_group.resseq: continue
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        yield_residue_group = False
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if atom_group.resname!=target_atom_group.resname: continue
          if verbose: print('    atom_group: altloc="%s" resname="%s"' % (
            atom_group.altloc, atom_group.resname))
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
        print(atom.name, end=' ')
      print(get_distance2(angle[0], angle[1]), end=' ')
      print(get_distance2(angle[1], angle[2]))
  return tmp

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
    print("tmp_key %s" % tmp_key)
    print("long_tmp_key %s" % long_tmp_key)
    print(atom1.quote())
    print(atom2.quote())
    print(is_n_glyco_bond(atom1, atom2))
    print(is_o_glyco_bond(atom1, atom2))
    print(is_glyco_bond(atom1, atom2))

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
        print('get_hand')
        print(c_atom, o_atom, angles)
      volume = get_chiral_volume(c_atom, o_atom, angles, verbose=verbose)
      hand = get_hand(c_atom, o_atom, angles, verbose=verbose) #"ALPHA"
      if hand is None:
        key = long_tmp_key
      else:
        data_link_key = "%s%s-%s" % (hand,
                                     c_atom.name.strip()[-1],
                                     o_atom.name.strip()[-1],
                                     )
        key = data_link_key
    else:
      print(" %s" % ("!"*86))
      _write_warning_line("  Possible link ignored")
      _write_warning_line(atom1.format_atom_record())
      _write_warning_line(atom2.format_atom_record())
      _write_warning_line("  N-linked glycan : %s" % (is_n_glyco_bond(atom1, atom2)))
      _write_warning_line("  O-linked glycan : %s" % (is_o_glyco_bond(atom1, atom2)))
      _write_warning_line("  Glycan-glycan   : %s" % (is_glyco_bond(atom1, atom2)))
      if c_atom is None: _write_warning_line("  No carbon atom found")
      if o_atom is None: _write_warning_line("  No oxygen atom found")
      print(" %s" % ("!"*86))
      #print " Distance", get_distance2(atom1, atom2)
      #assert 0
      #raise Sorry("Check input geometry")
      return None
  else:
    if verbose: print('bypass',long_tmp_key)
    key = long_tmp_key

  pdbres_pair = []
  for atom in [atom1, atom2]:
    pdbres_pair.append(atom.id_str(pdbres=True))
  if verbose:
    print("key2 %s" % key)
    print(pdbres_pair)
    print(atom1.quote())
    print(atom2.quote())
  return [pdbres_pair], [key], [(atom1, atom2)]

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

def check_for_acid(hierarchy, carbon):
  oxygens=0
  carbons=0
  bonds = get_bonded(hierarchy, carbon, bond_cutoff=1.8)
  for ba in bonds:
    if ba.element.strip()=="O":
      if len(get_bonded(hierarchy, ba, bond_cutoff=1.8))==1:
        oxygens+=1
  if carbon.element.strip()=="C":
    carbons=1
  if oxygens==2 and carbons==1:
    # acid
    return True
  return False

def check_valence_and_acid(hierarchy, atom):
  acid=None
  if atom.element.strip()=="C":
    acid = check_for_acid(hierarchy, atom)
  if acid: return False
  if atom.element.strip() not in linking_setup.simple_valence: return True
  bonds = get_bonded(hierarchy, atom, bond_cutoff=1.8)
  if len(bonds)==linking_setup.simple_valence[atom.element.strip()]:
    return False
  else:
    # check for acid
    if atom.element.strip()=="O":
      if len(bonds)==1:
        acid = check_for_acid(hierarchy, bonds[0])
        if acid: return False
  return True

def check_valence(hierarchy, atom):
  if atom.element.strip() not in linking_setup.simple_valence: return True
  bonds = get_bonded(hierarchy, atom, bond_cutoff=1.8)
  if len(bonds)==linking_setup.simple_valence[atom.element.strip()]:
    return False
  return True
