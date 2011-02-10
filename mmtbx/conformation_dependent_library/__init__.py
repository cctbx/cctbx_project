import sys
import copy
from string import letters, digits

import iotbx.pdb

from mmtbx.conformation_dependent_library.cdl_database import cdl_database

chararcters_36 = letters[:26]+digits

not_before_pro_groups = {
  "NonPGIV_nonxpro" : ["ALA",
                       "ARG",
                       "ASN",
                       "ASP",
                       "CYS",
                       "GLN",
                       "GLU",
                       "HIS",
                       "LEU",
                       "LYS",
                       "MET",
                       "PHE",
                       "SER",
                       "THR",
                       "TRP",
                       "TYR",
                       ],
  "IleVal_nonxpro" : ["ILE",
                      "VAL",
                      ],
  "Gly_nonxpro" : ["GLY"],
  "Pro_nonxpro" : ["PRO"],
}
before_pro_groups = {
  "NonPGIV_xpro" : not_before_pro_groups["NonPGIV_nonxpro"],
  "IleVal_xpro"  : not_before_pro_groups["IleVal_nonxpro"],
  "Gly_xpro"     : not_before_pro_groups["Gly_nonxpro"],
  "Pro_xpro"     : not_before_pro_groups["Pro_nonxpro"],
}
columns = [
  "",
  "",
  "mCNA", # C(-1) - N(0)  - Ca(0)
  "sCNA",
  "mNAB", # NAB   N(0)  - Ca(0) - Cb(0)
  "sNAB",
  "mNAC", # NAC   N(0)  - Ca(0) - C(0)
  "sNAC",
  "mBAC", # BAC   Cb(0) - Ca(0) - C(0)
  "sBAC",
  "mACO", # ACO   Ca(0) - C(0)  - O(0)
  "sACO",
  "mACN", # ACN   Ca(0) - C(0)  - N(+1)
  "sACN",
  "mOCN", # OCN   O(0)  - C(0)  - N(+1)
  "sOCN",
  "mCN",  # CN    C(-1) - N(0)
  "sCN",
  "mNA",  # NA    N(0)  - Ca(0)
  "sNA",
  "mAB",  # AB    Ca(0) - Cb(0)
  "sAB",
  "mAC",  # AC    Ca(0) - C(0)
  "sAC",
  "mCO",  # CO    C(0)  - O(0)
  "sCO",
  ]

def get_res_type_group(resname1, resname2):
  if resname2=="PRO":
    lookup = before_pro_groups
  else:
    lookup = not_before_pro_groups
  for key in lookup:
    if resname1 in lookup[key]:
      return key
  return None

def get_c_ca_n(atom_group):
  tmp = []
  for name in ["C", "CA", "N"]:
    for atom in atom_group.atoms():
      if atom.name.strip()==name:
        tmp.append(atom)
        break
    else:
      assert 0
  return tmp

def get_dihedral(xyzs):
  from elbow.chemistry.xyzClass import xyzClass
  for i, xyz in enumerate(xyzs):
    xyzs[i] = xyzClass(xyz)
  d = xyzs[0].BondDihedral(*xyzs[1:])
  return d

def round_to_ten(d):
  return int(round((float(d))/10))*10

def get_cdl_key(atom_group,
                last_atom_group,
                verbose=False,
                ):
  backbone = get_c_ca_n(atom_group)
  last_backbone = get_c_ca_n(last_atom_group)
  assert len(backbone)==3
  assert len(last_backbone)==3
  phi_atoms = [
    last_backbone[0],
    backbone[2],
    backbone[1],
    backbone[0],
    ]
  #phi_atoms.reverse()
  tmp = []
  for atom in phi_atoms:
    tmp.append(atom.xyz)
  phi = get_dihedral(tmp)
  psi_atoms = [
    last_backbone[2],
    last_backbone[1],
    last_backbone[0],
    backbone[2],
    ]
  #psi_atoms.reverse()
  tmp = []
  for atom in psi_atoms:
    tmp.append(atom.xyz)
  psi = get_dihedral(tmp)
  if verbose:
    print "psi, phi",psi,phi
  key = (round_to_ten(phi), round_to_ten(psi))
  return key

def get_i_seqs(threes):
  atoms = {}
  # i-1
  for name in ["C"]:
    for atom in threes[0].atoms():
      if atom.name.strip()==name:
        atoms["%s_minus_1" % name] = atom
        break
    else:
      assert 0
  # i
  for name in ["N", "CA", "CB", "C", "O"]:
    for atom in threes[1].atoms():
      if atom.name.strip()==name:
        atoms["%s_i" % name] = atom
        break
    else:
      if name!="CB":
        assert 0
  # i+1
  for name in ["N"]:
    for atom in threes[2].atoms():
      if atom.name.strip()==name:
        atoms["%s_plus_1" % name] = atom
        break
    else:
      assert 0
  return atoms

def apply_updates(threes,
                  restraint_values,
                  restraints_manager=None,
                  ):
  atoms = get_i_seqs(threes)
  for name in atoms:
    atom = atoms[name]
  for i, value in enumerate(restraint_values):
    if i<2: continue
    if columns[i][0]=="s": continue
    code = columns[i][1:]
    names = []
    if code=="CNA":   names = ["C_minus_1", "N_i",  "CA_i"      ]
    elif code=="NAB": names = ["N_i",       "CA_i", "CB_i"      ]
    elif code=="NAC": names = ["N_i",       "CA_i", "C_i"       ]
    elif code=="BAC": names = ["CB_i",      "CA_i", "C_i"       ]
    elif code=="ACO": names = ["CA_i",      "C_i",  "O_i"       ]
    elif code=="ACN": names = ["CA_i",      "C_i",  "N_plus_1"  ]
    elif code=="OCN": names = ["O_i",       "C_i",  "N_plus_1"  ]
    elif code=="CN":  names = ["C_i",  "N_plus_1" ]
    elif code=="NA":  names = ["C_i",  "CA_i" ]
    elif code=="AB":  names = ["CA_i", "CB_i" ]
    elif code=="AC":  names = ["CA_i", "C_i" ]
    elif code=="CO":  names = ["C_i",  "O_i" ]
    for j in range(len(names)):
      names[j] = atoms[names[j]].i_seq
    rnames = copy.deepcopy(names)
    rnames.reverse()
    if len(names)==3:
      for angle in restraints_manager.geometry.angle_proxies:
        if list(angle.i_seqs)==names or list(angle.i_seqs)==rnames:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            angle.i_seqs,
            angle.angle_ideal,
            angle.weight,
            restraint_values[i],
            restraint_values[i+1],
            )
          angle.angle_ideal = restraint_values[i]
          angle.weight = restraint_values[i+1]
          break
      else:
        print angle.i_seqs
        assert 0
    elif len(names)==2:
      bpt = restraints_manager.geometry.bond_params_table
      bond=None
      try:
        bond = bpt[names[0]][names[1]]
      except:
        try:
          bond = bpt[names[1]][names[0]]
        except:
          pass
      assert bond
      print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
        names,
        bond.distance_ideal,
        bond.weight,
        restraint_values[i],
        restraint_values[i+1],
        )
      bond.distance_ideal = restraint_values[i]
      bond.weight = restraint_values[i+1]
    else:
      assert 0

  assert 0

def update_restraints(hierarchy,
                      current_geometry=None, # xray_structure!!
                      restraints_manager=None,
                      dihedral_proxies=None,
                      verbose=False,
                      ):
  get_class = iotbx.pdb.common_residue_names_get_class
  #if verbose:
  #  hierarchy.show()
  threes = []
  for model in hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      for residue_group in chain.residue_groups():
        if verbose: print '  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode)

        threes.append(residue_group)
        while len(threes)>3:
          del threes[0]
        if len(threes)!=3: continue

        if residue_group.link_to_previous:
          for atom_group_i in threes[1].atom_groups():
            if get_class(atom_group_i.resname) not in ["common_amino_acid"]:
              continue
            if verbose:
              if atom_group_i.resname not in ["HOH"]:
                print '    atom_group_i: altloc="%s" resname="%s"' % (
                  atom_group_i.altloc, atom_group_i.resname)
            for atom_group_i_minus_1 in threes[0].atom_groups():
              if get_class(atom_group_i_minus_1.resname) not in ["common_amino_acid"]:
                continue
              if verbose: print '    atom_group_i_minus_1: altloc="%s" resname="%s"' % (
                atom_group_i_minus_1.altloc, atom_group_i_minus_1.resname)

              res_type_group = get_res_type_group(
                atom_group_i.resname,
                atom_group_i_minus_1.resname,
                )
              assert res_type_group
              print res_type_group, atom_group_i.resname, atom_group_i_minus_1.resname
              key = get_cdl_key(atom_group_i,
                                atom_group_i_minus_1,
                                verbose=verbose,
                                )
              tmp = (
                (key[0]+180)/10,
                (key[1]+180)/10,
                )
              #print tmp
              #print "%s%s" % (chararcters_36[tmp[0]-1],
              #                chararcters_36[tmp[1]-1])
              restraint_values = cdl_database[res_type_group][key]
              apply_updates(threes,
                            restraint_values,
                            restraints_manager,
                            )
              assert 0
          #assert 0
        else:
          for atom_group in residue_group.atom_groups():
            if verbose: print '    atom_group: altloc="%s" resname="%s"' % (
              atom_group.altloc, atom_group.resname)

  assert 0

#RESIDUE 20 ALA 21 PEPTIDE
#RESIDUE 21 LYS 22 PEPTIDE
#RESIDUE 22 GLU 23 PEPTIDE
test_string = """
RESIDUE 20 ALA 21 YAm6 21 XAjr
RESIDUE 21 LYS 22 YAjr 22 XEe7
RESIDUE 22 GLU 23 YEe7 23 XHl7
"""

def tst():
  for line in test_string.split("\n"):
    print line
    tmp = line.split()
    print tmp
  #assert 0

def run(filename):
  if False:
    for i in range(-188,188):
      print i,round_to_ten(i),abs(i-round_to_ten(i))
    assert 0
  if filename=="pdb3lo8.ent" and 0:
    tst()

  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  update_restraints(hierarchy,
                    verbose=True,
                    )

if __name__=="__main__":
  run(sys.argv[1])
