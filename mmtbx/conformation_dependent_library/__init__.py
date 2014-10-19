from __future__ import division
import sys
import copy
from string import letters, digits

from libtbx.utils import Sorry
import iotbx.pdb
from scitbx.math import dihedral_angle

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
headers = [
  "statistical type",
  "number",
  "C(-1) - N(0)  - Ca(0)",
  "",
  "N(0)  - Ca(0) - Cb(0)",
  "",
  "N(0)  - Ca(0) - C(0)",
  "",
  "Cb(0) - Ca(0) - C(0)",
  "",
  "Ca(0) - C(0)  - O(0)",
  "",
  "Ca(0) - C(0)  - N(+1)",
  "",
  "O(0)  - C(0)  - N(+1)",
  "",
  "C(-1) - N(0)",
  "",
  "N(0)  - Ca(0)",
  "",
  "Ca(0) - Cb(0)",
  "",
  "Ca(0) - C(0)",
  "",
  "C(0)  - O(0)",
  "",
  ]

def distance2(a,b):
  d2 = 0
  for i in range(3):
    d2 += (a.xyz[i]-b.xyz[i])**2
  return d2

def restraints_show(restraints_values):
  outl = ""
  for i, item in enumerate(restraints_values):
    if i%2==0:
      if i==0:
        s = "  %s, %s : %s %s\n"
      elif i<15:
        s = "  %-25s%s: %9.2f %9.2f\n"
      else:
        s = "  %-25s%s:   %9.4f %9.4f\n"
      outl += s % (headers[i],
                   headers[i+1],
                   restraints_values[i],
                   restraints_values[i+1],
        )
  return outl

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
  def __init__(self, geometry): #restraints_manager):
    self.geometry = geometry
    #self.restraints_manager = restraints_manager
    if geometry is None:
      self.bond_params_table = None
    else:
      self.bond_params_table = geometry.bond_params_table
    #except: self.bond_params_table = restraints_manager.bond_params_table
    self.errors = []

  def __repr__(self):
    return self.show()

  def show(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += " %s(%s)" % (residue.resname, residue.resseq)
    outl += " %s" % self.are_linked(return_value=True)
    return outl

  def show_detailed(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += "\nREMARK"
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def cis_group(self, limit=45., verbose=False):
    cis_peptide_bond = False
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue)
      ccn2, outl2 = get_c_ca_n(self[i-1])
      ca1 = ccn1[1]
      n = ccn1[2]
      c = ccn2[0]
      ca2 = ccn2[1]
      omega_atoms = [ca1, n, c, ca2]
      omega = dihedral_angle(sites=[atom.xyz for atom in omega_atoms], deg=True)
      if (180.-abs(omega))>limit:
        cis_peptide_bond = True
        break
    if verbose:
      if cis_peptide_bond:
        print 'cis peptide bond', cis_peptide_bond, omega
        print self
    return cis_peptide_bond

  def are_linked(self, return_value=False):
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue)
      ccn2, outl2 = get_c_ca_n(self[i-1])
      if ccn1 is None:
        for line in outl1:
          if line not in self.errors:
            self.errors.append(line)
        break
      if ccn2 is None:
        for line in outl2:
          if line not in self.errors:
            self.errors.append(line)
        break
      n = ccn1[2]
      c = ccn2[0]
      if self.bond_params_table is None:
        d2 = distance2(n,c)
        if d2<4: bond=True
        else: bond=False
      else:
        bond=self.bond_params_table.lookup(c.i_seq, n.i_seq)
      if not bond:
        #assert c.i_seq
        #assert n.i_seq
        break
    else:
      return True
    if return_value: return d2
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
    # i
    for name in ["N", "CA", "CB", "C", "O"]:
      for atom in self[1].atoms():
        if atom.name.strip()==name:
          atoms["%s_i" % name] = atom
          break
      #else:
      #  if name not in ["CB", "O"]:
      #    print self
      #    for atom in self[1].atoms():
      #      print atom.name, atom.xyz
      #    assert 0
    # i+1
    for name in ["N"]:
      for atom in self[2].atoms():
        if atom.name.strip()==name:
          atoms["%s_plus_1" % name] = atom
          break
    return atoms

  def get_cdl_key(self, exact=False, verbose=False):
    backbone_i_minus_1, junk = get_c_ca_n(self[0])
    backbone_i, junk = get_c_ca_n(self[1])
    backbone_i_plus_1, junk = get_c_ca_n(self[2])
    assert len(backbone_i_minus_1)==3
    assert len(backbone_i)==3
    assert len(backbone_i_plus_1)==3
    phi_atoms = [
      backbone_i_minus_1[0],
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      ]
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
    if exact: return (phi, psi)
    key = (round_to_ten(phi), round_to_ten(psi))
    return key

  def apply_updates(self,
                    restraint_values,
                    cdl_proxies,
                    ideal=True,
                    esd=True,
                    esd_factor=1.,
                    average=True,
                    verbose=False,
                    ):
    if not average:
      if restraint_values[0]=="I":
        print restraint_values
        assert 0
        return
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
        if angle_proxy is None:
          outl=""
          for key in atoms:
            outl += "\n    %-10s %s" % ( key, atoms[key].quote())
          raise Sorry("""CDL angle to be changed not set in model.
  Possible problems:
    Residue on special positions.

  Check:%s""" % outl)
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            angle_proxy.i_seqs,
            angle_proxy.angle_ideal,
            angle_proxy.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            )
        names.sort()
        registry[tuple(names)] = restraint_values[i]
        if ideal: angle_proxy.angle_ideal = restraint_values[i]
        if esd: angle_proxy.weight = esd_factor * 1/restraint_values[i+1]**2
      elif len(names)==2:
        bond=self.bond_params_table.lookup(*names)
        assert bond
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            names,
            bond.distance_ideal,
            bond.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            )
        names.sort()
        registry[tuple(names)] = restraint_values[i]
        #print "BOND", 1/restraint_values[i+1]**2/bond.weight,1/restraint_values[i+1]**2, bond.weight
        if ideal: bond.distance_ideal = restraint_values[i]
        if esd: bond.weight = esd_factor * 1/restraint_values[i+1]**2
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
#        for angle in self.restraints_manager.geometry.angle_proxies:
        for angle in self.geometry.angle_proxies:
          # could be better!
          akey = list(copy.deepcopy(angle.i_seqs))
          akey.sort()
          akey=tuple(akey)
          if akey==key or akey==rkey:
            angle.angle_ideal = averages[key]/averages.n[key]
            break
        else:
          print key,rkey
          print averages[key]
          print averages.n[key]
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
  outl = []
  for name in ["C", "CA", "N"]:
    for atom in atom_group.atoms():
      if atom.name.strip()==name:
        tmp.append(atom)
        break
    else:
      for atom in atom_group.atoms():
        outl.append(atom.format_atom_record())
      tmp = None
      break
  return tmp, outl

def round_to_ten(d):
  t = int(round((float(d))/10))*10
  if t==180: return -180
  return t

def get_restraint_values(threes, interpolate=False):
  from mmtbx.conformation_dependent_library import utils
  res_type_group = get_res_type_group(
    threes[1].resname,
    threes[2].resname,
  )
  if res_type_group is None: return None
  if interpolate:
    restraint_values = ["2", -1]
    key = threes.get_cdl_key(exact=interpolate)
    for i in range(2,26):
      grid = utils.get_grid_values(res_type_group, key[0], key[1], column=i)
      #utils.print_grid(grid, key[0], key[1])
      index = utils.get_index(*key)
      r = utils.interpolate_2d(grid, index)
      #print i, key, index, r
      restraint_values.append(r)
  else:
    key = threes.get_cdl_key()
    restraint_values = cdl_database[res_type_group][key]
  return restraint_values

def generate_protein_threes(hierarchy,
                            geometry, #restraints_manager,
                            include_non_linked=False,
                            verbose=False,
                            ):
  get_class = iotbx.pdb.common_residue_names_get_class
  threes = ThreeProteinResidues(geometry) #restraints_manager)
  for model in hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      for conformer in chain.conformers():
        if verbose: print '  conformer: altloc="%s"' % (
          conformer.altloc)
        while threes: del threes[0]
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print '    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid())
              #for atom in residue.atoms():
              #  if verbose: print '         atom: name="%s"' % (atom.name)
          if verbose:
            print 'residue class : %s' % get_class(residue.resname)
          if get_class(residue.resname) not in ["common_amino_acid"]:
            continue
          if include_non_linked:
            list.append(threes, residue)
            if len(threes)>3: del threes[0]
          else:
            threes.append(residue)
          if len(threes)!=3: continue
          assert len(threes)<=3
          yield threes
      threes = ThreeProteinResidues(geometry)

def setup_restraints(geometry, # restraints_manager
                     verbose=False,
                     ):
  ba_registry = bond_angle_registry()
  for angle in geometry.angle_proxies:
    ba_registry[angle.i_seqs]=angle
  return ba_registry

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      cdl_proxies=None,
                      ideal=True,
                      esd=True,
                      esd_factor=1.,
                      interpolate=False,
                      log=None,
                      verbose=False,
                      ):
  global registry
  registry = RestraintsRegistry()
  if current_geometry:
    assert not sites_cart
    sites_cart = current_geometry.sites_cart()
  if sites_cart:
    pdb_atoms = hierarchy.atoms()
    #if atom_lookup:
    #  for j_seq, scatterer in enumerate(current_geometry.scatterers()):
    #    pdb_atoms[atom_lookup[scatterer.label]].xyz = sites_cart[j_seq]
    #else:
    # XXX PDB_TRANSITION VERY SLOW
    for j_seq, atom in enumerate(pdb_atoms):
      atom.xyz = sites_cart[j_seq]
      #atom_lookup[atom.id_str()] = j_seq

  threes = None
  average_updates = 0
  total_updates = 0
  for threes in generate_protein_threes(hierarchy,
                                        geometry, #restraints_manager,
                                        #verbose=verbose,
                                        ):
    if threes.cis_group():
      if verbose and 0:
        print 'cis '*20
        print threes
      continue

    if 0:
      res_type_group = get_res_type_group(
        threes[1].resname,
        threes[2].resname,
         )
      if res_type_group is None: continue
      key = threes.get_cdl_key() #verbose=verbose)
      restraint_values = cdl_database[res_type_group][key]
      print restraint_values
      print len(restraint_values)
      assert 0
    else:
      restraint_values = get_restraint_values(threes, interpolate=interpolate)

    #if 1:
    #  print threes, threes.are_linked(), res_type_group, key, restraint_values

    if restraint_values is None: continue

    if restraint_values[0]=="I":
      #print threes, threes.are_linked(), res_type_group, key, restraint_values[:4]
      average_updates += 1
    else:
      total_updates += 1
    threes.apply_updates(restraint_values,
                         cdl_proxies,
                         ideal=ideal,
                         esd=esd,
                         esd_factor=esd_factor,
                         )
  if registry.n: threes.apply_average_updates(registry)
#  restraints_manager.geometry.reset_internals()
  geometry.reset_internals()
  if verbose and threes and threes.errors:
    if log:
      log.write("  Residues not completely updated with CDL restraints\n\n")
    for line in threes.errors:
      if log:
        log.write("%s\n" % line)
      else:
        print line
#  print 'average updates',average_updates,total_updates
#  assert average_updates==0
  return geometry #restraints_manager

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
  if 0:
    psi = -180
    lookup = "Gly_nonxpro"
    print lookup
    for phi in range(170,181):
      key = (round_to_ten(psi),round_to_ten(phi))
      print 'key',psi,phi,round_to_ten(psi),round_to_ten(phi),key,
      print cdl_database[lookup][key][:4]

  run(sys.argv[1])
