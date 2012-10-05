import iotbx

from scitbx.array_family import flex
from mmtbx.chemical_components import get_type

get_class = iotbx.pdb.common_residue_names_get_class

sugar_types = ["SACCHARIDE",
               "D-SACCHARIDE",
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

class empty:
  def __repr__(self):
    outl = ""
    for attr in self.__dict__:
      outl += "  %s : %s\n" % (attr, getattr(self, attr))
    return outl

def get_distance2(atom1, atom2):
  d2 = (atom1.xyz[0]-atom2.xyz[0])**2
  d2 += (atom1.xyz[1]-atom2.xyz[1])**2
  d2 += (atom1.xyz[2]-atom2.xyz[2])**2
  return d2

def get_chiral_volume(centre, atom1, atom2, atom3):
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

def is_glyco_bond(atom1, atom2):
  if get_type(atom1.parent().resname) is None: return False
  if get_type(atom2.parent().resname) is None: return False
  if not get_type(atom1.parent().resname).upper() in sugar_types: return False
  if not get_type(atom2.parent().resname).upper() in sugar_types: return False
  return True

def is_glyco_amino_bond(atom1, atom2):
  if get_type(atom1.parent().resname) is None: return False
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

def get_hand(c_atom, o_atom, angles):
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
    print '-'*80
    for atom in others:
      print atom.format_atom_record()
    return None
  v = get_chiral_volume(*others)
  if v<0:
    return "BETA"
  else:
    return "ALPHA"

def get_classes(atom, verbose=False):
  attrs = [
    "common_sugar", # not in get_class
    "common_water",
    "common_element",
    "common_small_molecule",
    "common_amino_acid",
    "common_rna_dna",
    "other",
    "unknown",
    ]
  atom_group = atom.parent()
  classes = empty()
  for attr in attrs:
    setattr(classes, attr, False)
  if verbose:
    print '    atom_group1: altloc="%s" resname="%s" class="%s"' % (
      atom_group.altloc,
      atom_group.resname,
      get_class(atom_group.resname),
      )
  gc = get_class(atom_group.resname)
  for i, attr in enumerate(attrs):
    rc = None
    if i:
      rc = gc
    else:
      if(get_type(atom_group.resname) is not None and
         get_type(atom_group.resname).upper() in sugar_types):
        rc = attr
    if rc==attr:
      setattr(classes, attr, True)
  return classes

def get_closest_atoms(atom_group1,
                      atom_group2,
                      ignore_hydrogens=True,
                      ):
  min_d2 = 1e5
  min_atom1 = None
  min_atom2 = None
  for i, atom1 in enumerate(atom_group1.atoms()):
    if ignore_hydrogens:
      if atom1.element.strip() in ["H", "D"]: continue
    for j, atom2 in enumerate(atom_group2.atoms()):
      if ignore_hydrogens:
        if atom2.element.strip() in ["H", "D"]: continue
      #if i>=j: continue
      d2 = get_distance2(atom1, atom2)
      if d2<min_d2:
        min_atom1 = atom1
        min_atom2 = atom2
        min_d2 = d2
  return min_atom1, min_atom2

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

def process_atom_groups_for_linking(pdb_hierarchy,
                                    atom1,
                                    atom2,
                                    classes1,
                                    classes2,
                                    bond_cutoff=2.,
                                    verbose=False,
                                    ):
  bond_cutoff *= bond_cutoff
  atom_group1 = atom1.parent()
  atom_group2 = atom2.parent()
  residue_group1 = atom_group1.parent()
  residue_group2 = atom_group2.parent()
  atom1, atom2 = get_closest_atoms(residue_group1, residue_group2)
  if get_distance2(atom1, atom2)>bond_cutoff: return None
  if is_glyco_bond(atom1, atom2):
    # glyco bonds need to be in certain order
    if atom1.name.find("C")>-1:
      tmp_atom = atom1
      atom1 = atom2
      atom2 = tmp_atom
    
  elif is_glyco_amino_bond(atom1, atom2):
    # problem in 3sgk 
#------------------------------------------------------------------
#ATOM   3803  OD1 ASP C  64      19.148  52.821 -19.425  1.00 70.10
#HETATM 5030  O6  NAG C1461      19.450  52.248 -18.258  1.00 60.78
#distance 1.33469921705
#------------------------------------------------------------------
    if atom1.element.strip()=="O" and atom2.element.strip()=="O": return None
    if atom2.name.find("C")>-1: # needs to be better uding get_class???
      tmp_atom = atom1
      atom1 = atom2
      atom2 = tmp_atom

  tmp_key = "%s-%s" % (atom1.parent().resname.strip(),
                       atom2.parent().resname.strip(),
                       )

  if is_n_glyco_bond(atom1, atom2):
    if tmp_key in standard_n_links:
      data_links = ""
    key = tmp_key
  elif is_o_glyco_bond(atom1, atom2):
    if tmp_key in standard_o_links:
      data_links = ""
    key = tmp_key
  elif is_glyco_bond(atom1, atom2):
    data_links = ""
    c_atom = None
    o_atom = None
    if atom1.name.find("C")>-1:
      c_atom = atom1
    elif atom2.name.find("C")>-1:
      c_atom = atom2
    if atom1.name.find("O")>-1:
      o_atom = atom1
    elif atom2.name.find("O")>-1:
      o_atom = atom2
    if c_atom and o_atom:
      angles = get_angles_from_included_bonds(pdb_hierarchy,
                                              [[atom1, atom2]],
                                              bond_cutoff=2., #control.intra_residue_bond_cutoff,
                                              )
      hand = get_hand(c_atom, o_atom, angles) #"ALPHA"
      assert hand
        
      data_link_key = "%s%s-%s" % (hand,
                                   c_atom.name.strip()[-1],
                                   o_atom.name.strip()[-1],
                                   )
      if data_link_key in [
        "BETAB-B",
        ]: assert 0
      #cif_links = cif_links.replace(tmp_key, data_link_key)
      key = data_link_key
    else:
      print " %s" % ("!"*84)
      print _write_warning_line("  Possible link ignored")
      print _write_warning_line(atom1.format_atom_record())
      print _write_warning_line(atom2.format_atom_record())
      print _write_warning_line("  N-linked glycan : %s" % (is_n_glyco_bond(atom1, atom2)))
      print _write_warning_line("  O-linked glycan : %s" % (is_o_glyco_bond(atom1, atom2)))
      print _write_warning_line("  Glycan-glycan   : %s" % (is_glyco_bond(atom1, atom2)))
      if c_atom is None: print _write_warning_line("  No carbon atom found")
      if o_atom is None: print _write_warning_line("  No oxygen atom found")
      print " %s" % ("!"*84)
  else:
    key = tmp_key

  pdbres_pair = []
  for atom in [atom1, atom2]:
    pdbres_pair.append(atom.id_str(pdbres=True))
  return pdbres_pair, key, (atom1, atom2)
