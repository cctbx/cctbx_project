from __future__ import absolute_import, division, print_function
import scitbx.math
from libtbx.utils import Sorry
import libtbx.load_env
from libtbx import group_args
import iotbx.pdb
from cctbx.array_family import flex
import mmtbx.rotamer
import sys, os
from six.moves import zip
from six.moves import range

def find_source_dir(optional=False):
  result = libtbx.env.find_in_repositories(os.path.join("mmtbx", "rotamer"))
  if result is None and not optional:
    raise Sorry("""\
Can't seem to find mmtbx/rotamer/ directory.
  """)
  return result

class PropertyFile:

  #properties = {}

  def __init__(self):
    self.properties = {}

  def process(self, fileLoc):
    try:
      f = open(fileLoc)
    except ImportError as e:
      print(fileLoc+" file not found")
      sys.exit()
    for line in f:
      if (line.startswith("#") or line == "\n"): continue
      else:
        props = line.split("=")
        self.properties[props[0].strip()] = props[1].strip().strip("\"")
    f.close()

class SidechainAngles:

  #knownAA = {}
  chisPerAA = {}
  anglesForAA = {}
  atomsForAngle = {}
  rotamersForAA = {}
  anglesForRot = {}
  atomsMoveWithAngle = {}
  resAtomsToChi = {}
  frequencies_from_rotamer = {}

  def __init__(self, show_errs):
    self.show_errors = show_errs
    source_dir = find_source_dir()
    #print source_dir
    f = PropertyFile()
    f.process(os.path.join(source_dir, "sidechain_angles.props"))
    for aa in f.properties['aminoacids'].split(","):
      #print aa + f.properties[aa+".chis"], f.properties['%s.frequencies' % aa]
      for rot, freq in zip(f.properties['%s.rotamers' % aa].split(','),
                           f.properties['%s.frequencies' % aa].split(','),
                           ):
        if not freq: continue
        self.frequencies_from_rotamer.setdefault(aa, {})
        self.frequencies_from_rotamer[aa][rot] = freq
      self.resAtomsToChi[aa] = {}
      self.chisPerAA[aa] = f.properties[aa+".chis"] #gives aaName -> # of chis
      #print f.properties[aa+".angles"].split(",")
      anglelist = f.properties[aa+".angles"].split(",")
      rotamerlist = f.properties[aa+".rotamers"].split(",")
      #print anglelist.count('')
      self.anglesForAA[aa] = anglelist #aaName -> [mobile angles]
      chi_ctr = 0
      for angle in anglelist:
        if angle != '':
          key = aa+"."+angle
          #print key
          self.atomsForAngle[key] = f.properties[key].split(",") #aaName.angle -> atoms
          chi_atoms = f.properties[key]
          if aa == 'val' and angle == 'chi1':
            chi_atoms = chi_atoms.replace('CG1', 'CG2')
          elif aa == 'thr' and angle == 'chi1':
            chi_atoms = chi_atoms.replace('OG1', 'CG2')
          elif aa == 'ile' and angle == 'chi1':
            chi_atoms = chi_atoms.replace('CG1', 'CG2')
          self.resAtomsToChi[aa][chi_atoms] = angle
          if aa == 'leu' or aa == 'val' or aa == 'thr' or aa == 'arg':
            if chi_ctr < int(f.properties[aa+".chis"]):
              key2 = key + "_atoms"
              #print f.properties[aa+".chis"]
              #print chi_ctr
              self.atomsMoveWithAngle[key] = f.properties[key2].split(",")
        chi_ctr+=1
      self.rotamersForAA[aa] = rotamerlist
      for rotamer in rotamerlist:
        if rotamer != '':
          key = aa+"."+rotamer
          self.anglesForRot[key] = f.properties[key].split(" ")

  def get_rotamers(self, residue_name):
    rotamers = {}
    aa = residue_name.lower()
    if (not aa in self.rotamersForAA):
      return None
    rotamer_list = self.rotamersForAA.get(aa)
    for rotamer in rotamer_list :
      if (rotamer != ""):
        key = aa + "." + rotamer
        rotamers[rotamer] = [ float(x) for x in self.anglesForRot[key] ]
    return rotamers

  def get_rotamer_angles(self, residue_name, rotamer_name):
    return_angles = []
    aa = residue_name.lower()
    key = aa + '.' + rotamer_name
    angles = self.anglesForRot.get(key)
    if angles is None:
      return None
    for angle in angles:
      return_angles.append(float(angle))
    return return_angles

  def get_rotamer_expectation_frequencies(self, residue_name, rotamer_name):
    aa = residue_name.lower()
    if (not aa in self.rotamersForAA): return None
    freqs = self.frequencies_from_rotamer.get(aa)
    if not freqs: return None
    return freqs.get(rotamer_name, None)

  def measureChiAngles(
        self,
        res,
        atom_dict = None,
        sites_cart = None):
    resName = res.resname.lower().strip()
    get_class = iotbx.pdb.common_residue_names_get_class
    if(get_class(res.resname) == "common_amino_acid"):
      try:
        numChis = int(self.chisPerAA[resName])
        values = []
        for i in range(numChis):
          values.append(self.measureAngle(
                          angleName="chi"+str(i+1),
                          res=res,
                          atom_dict=atom_dict,
                          sites_cart=sites_cart))
        return values
      except KeyError:
        if self.show_errors: print(resName + " is an unknown residue type")
        return None
    else:
      return None

#  def measureChiAnglesDict(self, atom_dict, resName):
#    try:
#      numChis = int(self.chisPerAA[resName])
#      values = []
#      for i in range(numChis):
#        values.append(self.measureAngle("chi"+str(i+1), res))
#      return values
#    except KeyError:
#      resName + " is unknown"

  def measureAngle(self, angleName, res, atom_dict, sites_cart=None):
    angleAtoms = self.extract_chi_atoms(
                   angleName=angleName,
                   res=res,
                   atom_dict=atom_dict)
    if angleAtoms is None:
      return None
    if sites_cart is None:
      sites = [a.xyz for a in angleAtoms]
    else:
      sites = []
      for a in angleAtoms:
        sites.append(sites_cart[a.i_seq])
    return scitbx.math.dihedral_angle(
      sites=sites, deg=True)

  def extract_chi_atoms(self, angleName, res, atom_dict=None):
    atomNamesMap = None
    if (atom_dict is None):
      atomNamesMap = makeAtomDict(res)
    else:
      atomNamesMap = atom_dict
    atomNames = self.atomsForAngle[res.resname.lower().strip()+"."+angleName]
    angleAtoms = []
    for at in atomNames:
      namelist = at.split(";")
      testAtom = None
      j = 0
      while (testAtom == None and j < len(namelist)):
        testAtom = atomNamesMap.get(namelist[j])
        j += 1
      if (testAtom != None):
        angleAtoms.append(atomNamesMap.get(testAtom.name))
      else:
        return None
    return angleAtoms

def makeAtomDict(res):
  atomNamesMap = {}
  for atom in res.atoms():
    atomNamesMap[atom.name] = atom
  return atomNamesMap

def collect_sidechain_chi_angles(pdb_hierarchy, atom_selection=None):

  angle_lookup = SidechainAngles(False)
  residue_chis = []
  if atom_selection is not None:
    if (isinstance(atom_selection, flex.bool)):
      actual_selection = atom_selection
    elif (isinstance(atom_selection, flex.size_t)):
      actual_selection = flex.bool(pdb_hierarchy.atoms_size(), False)
      actual_selection.set_selected(atom_selection, True)
  if atom_selection is None:
    actual_selection = flex.bool(pdb_hierarchy.atoms_size(), True)

  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          n_chi = angle_lookup.chisPerAA.get(residue.resname.lower(), 0)
          try :
            n_chi = int(n_chi)
          except ValueError :
            continue
          chis = []
          altloc = residue.atoms()[0].fetch_labels().altloc
          i_seqs = []
          for i in range(1, n_chi+1):
            atoms = angle_lookup.extract_chi_atoms("chi%d" % i, residue)
            if atoms is None:
              pass
            else :
              i_seqs = [ atom.i_seq for atom in atoms ]
              if actual_selection.select(flex.size_t(i_seqs)).all_eq(True):
                chis.append(group_args(chi_id=i, i_seqs=i_seqs))
          if len(chis) > 0:
            residue_info = group_args(
              residue_name=residue.resname,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              chis=chis)
            residue_chis.append(residue_info)
  return residue_chis

def collect_residue_torsion_angles(pdb_hierarchy,
                                    atom_selection=None,
                                    chi_angles_only=False):
  get_class = iotbx.pdb.common_residue_names_get_class
  residue_torsions = []

  ### chi angles ###
  residue_chis = collect_sidechain_chi_angles(
                   pdb_hierarchy=pdb_hierarchy,
                   atom_selection=atom_selection)
  residue_torsions = residue_chis
  if chi_angles_only:
    return residue_torsions

  ##################

  if atom_selection is not None:
    if (isinstance(atom_selection, flex.bool)):
      actual_selection = atom_selection
    elif (isinstance(atom_selection, flex.size_t)):
      actual_selection = flex.bool(pdb_hierarchy.atoms_size(), False)
      actual_selection.set_selected(atom_selection, True)
  if atom_selection is None:
    actual_selection = flex.bool(pdb_hierarchy.atoms_size(), True)
  previous_residue = None
  next_residue = None
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for i_res, residue in enumerate(conformer.residues()):
          if (get_class(residue.resname) != "common_amino_acid"):
            continue
          if i_res < (len(conformer.residues())-1):
            next_residue = conformer.residues()[i_res+1]
          else:
            next_residue = None
          torsions = []
          # atoms_to_work = [prevCA, prevC, curN, curCA, curC, nextN]
          atoms_to_work = [None]*6
          atoms_to_work[2] = residue.find_atom_by(name=" N  ")
          atoms_to_work[3] = residue.find_atom_by(name=" CA ")
          atoms_to_work[4] = residue.find_atom_by(name=" C  ")
          if previous_residue is not None:
            atoms_to_work[0] = previous_residue.find_atom_by(name=" CA ")
            atoms_to_work[1] = previous_residue.find_atom_by(name=" C  ")
          if next_residue is not None:
            atoms_to_work[5] = next_residue.find_atom_by(name=" N  ")

          # atoms_to_work = [prevCA, prevC, curN, curCA, curC, nextN]
          for i in range(len(atoms_to_work)):
            if (atoms_to_work[i] is not None and
                not actual_selection[atoms_to_work[i].i_seq]):
              atoms_to_work[i] = None
          for i, name in enumerate(["omega", "phi", "psi"]):
            if atoms_to_work[i:i+4].count(None) == 0:
              angle = mmtbx.rotamer.omega_from_atoms(
                  atoms_to_work[i],
                  atoms_to_work[i+1],
                  atoms_to_work[i+2],
                  atoms_to_work[i+3])
              if angle is not None:
                i_seqs = [
                    atoms_to_work[i].i_seq,
                    atoms_to_work[i+1].i_seq,
                    atoms_to_work[i+2].i_seq,
                    atoms_to_work[i+3].i_seq]
                torsions.append(group_args(chi_id=name, i_seqs=i_seqs))
          altloc = residue.atoms()[0].fetch_labels().altloc
          if len(torsions) > 0:
            residue_info = group_args(
              residue_name=residue.resname,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              chis=torsions)
            residue_torsions.append(residue_info)
          previous_residue = residue
  return residue_torsions
