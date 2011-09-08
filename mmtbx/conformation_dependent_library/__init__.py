import sys
import copy
from string import letters, digits

import iotbx.pdb

from mmtbx.conformation_dependent_library.cdl_database import cdl_database
from mmtbx.conformation_dependent_library.bond_angle_registry import \
  bond_angle_registry

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

class RestraintsRegistry(dict):
  def __init__(self):
    self.n = {}

  def __setitem__(self, key, item):
    if key in self:
      if self[key]!=item:
        self.n.setdefault(key,1)
        self.n[key]+=1
        dict.__setitem__(self, key, (self[key]+item))
    else:
      dict.__setitem__(self, key, item)

registry = RestraintsRegistry()

class ThreeProteinResidues(list):
  def __init__(self, restraints_manager):
    self.restraints_manager = restraints_manager
    self.bond_params_table = restraints_manager.geometry.bond_params_table
    #except: self.bond_params_table = restraints_manager.bond_params_table

  def __repr__(self):
    return self.show()

  def show(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += " %s(%s)" % (residue.resname, residue.resseq)
    return outl

  def show_detailed(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += "\nREMARK"
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def are_linked(self):
    for i, residue in enumerate(self):
      if i==0: continue
      try: ccn1 = get_c_ca_n(residue)
      except Exception: break
      try: ccn2 = get_c_ca_n(self[i-1])
      except Exception: break
      n = ccn1[2]
      c = ccn2[0]
      bond=self.bond_params_table.lookup(c.i_seq, n.i_seq)
      if not bond:
        #assert c.i_seq
        #assert n.i_seq
        break
    else:
      return True
    return False

  def append(self, residue):
    list.append(self, residue)
    while len(self)>3:
      del self[0]
    if len(self)>=2:
      while not self.are_linked():
        del self[0]
        if len(self)==0: break

  def get_i_seqs(self):
    atoms = {}
    # i-1
    for name in ["C"]:
      for atom in self[0].atoms():
        if atom.name.strip()==name:
          atoms["%s_minus_1" % name] = atom
          break
      else:
        assert 0
    # i
    for name in ["N", "CA", "CB", "C", "O"]:
      for atom in self[1].atoms():
        if atom.name.strip()==name:
          atoms["%s_i" % name] = atom
          break
      else:
        if name not in ["CB", "O"]:
          print self
          for atom in self[1].atoms():
            print atom.name, atom.xyz
          assert 0
    # i+1
    for name in ["N"]:
      for atom in self[2].atoms():
        if atom.name.strip()==name:
          atoms["%s_plus_1" % name] = atom
          break
      else:
        assert 0
    return atoms

  def get_cdl_key(self, verbose=False):
    backbone_i_minus_1 = get_c_ca_n(self[0])
    backbone_i = get_c_ca_n(self[1])
    backbone_i_plus_1 = get_c_ca_n(self[2])
    assert len(backbone_i_minus_1)==3
    assert len(backbone_i)==3
    assert len(backbone_i_plus_1)==3
    phi_atoms = [
      backbone_i_minus_1[0],
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      ]
    from scitbx.math import dihedral_angle
    phi = dihedral_angle(sites=[atom.xyz for atom in phi_atoms], deg=True)
    psi_atoms = [
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      backbone_i_plus_1[2],
      ]
    psi = dihedral_angle(sites=[atom.xyz for atom in psi_atoms], deg=True)
    if verbose:
      print "psi, phi",psi,phi
    key = (round_to_ten(phi), round_to_ten(psi))
    return key

  def apply_updates(self,
                    restraint_values,
                    cdl_proxies,
                    ideal=True,
                    esd=True,
                    verbose=False,
                    ):
    atoms = self.get_i_seqs()
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
      elif code=="CN":  names = ["C_minus_1",  "N_i" ]
      elif code=="NA":  names = ["N_i",  "CA_i" ]
      elif code=="AB":  names = ["CA_i", "CB_i" ]
      elif code=="AC":  names = ["CA_i", "C_i" ]
      elif code=="CO":  names = ["C_i",  "O_i" ]
      # not all amino acids have a CB
      if "CB_i" in names and not "CB_i" in atoms: continue
      # sometimes the O is not in the model
      if "O_i" in names and not "O_i" in atoms: continue
      for j in range(len(names)):
        names[j] = atoms[names[j]].i_seq
      if len(names)==3:
        rnames = copy.deepcopy(names)
        rnames.reverse()
        angle_proxy = cdl_proxies.get(tuple(names), None)
        if angle_proxy is None:
          angle_proxy = cdl_proxies.get(tuple(rnames), None)
        assert angle_proxy
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            angle_proxy.i_seqs,
            angle_proxy.angle_ideal,
            angle_proxy.weight,
            restraint_values[i],
            restraint_values[i+1],
            )
          #print "ANGLE", 1/restraint_values[i+1]**2/angle.weight,1/restraint_values[i+1]**2, angle.weight
        if ideal: angle_proxy.angle_ideal = restraint_values[i]
        if esd: angle_proxy.weight = 1/restraint_values[i+1]**2
      elif len(names)==2:
        bond=self.bond_params_table.lookup(*names)
        assert bond
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            names,
            bond.distance_ideal,
            bond.weight,
            restraint_values[i],
            restraint_values[i+1],
            )
        names.sort()
        registry[tuple(names)] = restraint_values[i]
        #print "BOND", 1/restraint_values[i+1]**2/bond.weight,1/restraint_values[i+1]**2, bond.weight
        if ideal: bond.distance_ideal = restraint_values[i]
        if esd: bond.weight = 1/restraint_values[i+1]**2
      else:
        assert 0

  def apply_average_updates(self, averages, verbose=False):
    if verbose:
      print averages
      print averages.n
    for key in averages.n:
      if len(key)==2:
        bond=self.bond_params_table.lookup(*key)
        bond.distance_ideal = averages[key]/averages.n[key]
      elif len(key)==3:
        rkey = list(copy.deepcopy(key))
        rkey.reverse()
        rkey=tuple(rkey)
        print averages
        print averages.n
        for angle in self.restraints_manager.geometry.angle_proxies:
          if angle.i_seqs==key or angle.i_seqs==rkey:
            print angle.angle_ideal,
            angle.angle_ideal = averages[key]
            print angle.angle_ideal
            assert 0
            break
        else:
          assert 0

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
      print "Residues not completely updated with CDL restraints"
      for atom in atom_group.atoms():
        print atom.format_atom_record()
      assert 0
  return tmp

def round_to_ten(d):
  t = int(round((float(d))/10))*10
  if t==180: return -180
  return t

def generate_protein_threes(hierarchy,
                            restraints_manager,
                            verbose=False,
                            ):
  get_class = iotbx.pdb.common_residue_names_get_class
  threes = ThreeProteinResidues(restraints_manager)
  for model in hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      for conformer in chain.conformers():
        if verbose: print '  conformer: altloc="%s" only_residue="%s"' % (
          conformer.altloc, conformer.only_residue)
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print '    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid())
              for atom in residue.atoms():
                if verbose: print '         atom: name="%s"' % (atom.name)
          if verbose:
            print 'residue class : %s' % get_class(residue.resname)
          if get_class(residue.resname) not in ["common_amino_acid"]:
            continue
          threes.append(residue)
          if len(threes)!=3: continue
          assert len(threes)<=3
          yield threes

def setup_restraints(restraints_manager,
                     verbose=False,
                     ):
  ba_registry = bond_angle_registry()
  for angle in restraints_manager.geometry.angle_proxies:
    ba_registry[angle.i_seqs]=angle
  return ba_registry

def update_restraints(hierarchy,
                      restraints_manager,
                      current_geometry=None, # xray_structure!!
                      cdl_proxies=None,
                      ideal=True,
                      esd=True,
                      verbose=False,
                      ):
  global registry
  registry = RestraintsRegistry()
  if current_geometry:
    sites_cart = current_geometry.sites_cart()
    pdb_atoms = hierarchy.atoms()
    # XXX PDB_TRANSITION SLOW
    for j_seq,atom in enumerate(pdb_atoms):
      atom.xyz = sites_cart[j_seq]

  for threes in generate_protein_threes(hierarchy,
                                        restraints_manager,
                                        verbose=verbose,
                                        ):
    res_type_group = get_res_type_group(
      threes[1].resname,
      threes[2].resname,
      )
    if res_type_group is None:
      #print "Non standard amino-acid skipped"
      #print threes
      continue

    key = threes.get_cdl_key(verbose=verbose)
    if verbose: print 'key',key

    restraint_values = cdl_database[res_type_group][key]
    if verbose:
      print threes, threes.are_linked(), res_type_group, key, restraint_values

    threes.apply_updates(restraint_values,
                         cdl_proxies,
                         ideal=ideal,
                         esd=esd,
                         )
  if registry.n: threes.apply_average_updates(registry)
  restraints_manager.geometry.reset_internals()
  return restraints_manager

def run(filename):
  if False:
    for i in range(-188,188):
      print i,round_to_ten(i),abs(i-round_to_ten(i))
    assert 0

  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  update_restraints(hierarchy,
                    #verbose=True,
                    )

if __name__=="__main__":
  run(sys.argv[1])
