import sys
from string import letters, digits

import iotbx.pdb

from elbow.scripts.cdl_database import cdl_database

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

def update_restraints(hierarchy,
                      current_geometry=None,
                      dihedral_proxies=None,
                      verbose=False,
                      ):
  get_class = iotbx.pdb.common_residue_names_get_class
  if verbose:
    hierarchy.show()
  for model in hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      last_residue_group = None
      for residue_group in chain.residue_groups():
        if verbose: print '  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode)
        if residue_group.link_to_previous:
          for atom_group in residue_group.atom_groups():
            if verbose:
              if atom_group.resname not in ["HOH"]:
                print '    atom_group: altloc="%s" resname="%s"' % (
                  atom_group.altloc, atom_group.resname)
            if atom_group.altloc.strip()!="": continue
            if get_class(atom_group.resname) in ["common_amino_acid"]:
              for last_atom_group in last_residue_group.atom_groups():
                if verbose: print '    last_atom_group: altloc="%s" resname="%s"' % (
                  last_atom_group.altloc, last_atom_group.resname)
                if last_atom_group.altloc.strip()!="": continue
                res_type_group = get_res_type_group(atom_group.resname,
                                                    last_atom_group.resname,
                                                    )
                print 'res_type_group',atom_group.resname,last_atom_group.resname,res_type_group
                assert res_type_group
                #print cdl_database[res_type_group]
                key = get_cdl_key(atom_group, last_atom_group, verbose=verbose)
                print key
                #print chararcters_36, letters, digits
                tmp = (
                  (key[0]+180)/10,
                  (key[1]+180)/10,
                  )
                print tmp
                print "%s%s" % (chararcters_36[tmp[0]-1],
                                chararcters_36[tmp[1]-1])
                #print cdl_database[res_type_group][key]
          #assert 0
        else:
          for atom_group in residue_group.atom_groups():
            if verbose: print '    atom_group: altloc="%s" resname="%s"' % (
              atom_group.altloc, atom_group.resname)

        last_residue_group = residue_group
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
  update_restraints(hierarchy,
                    verbose=True,
                    )

if __name__=="__main__":
  run(sys.argv[1])
