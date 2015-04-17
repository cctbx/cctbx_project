from __future__ import division
import scitbx.math
from libtbx.utils import Sorry
import libtbx.load_env
from libtbx import group_args
import iotbx.pdb
from cctbx.array_family import flex
import mmtbx.rotamer
import sys, os

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
    except ImportError, e:
      print fileLoc+" file not found"
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

  def __init__(self, show_errs):
    self.show_errors = show_errs
    source_dir = find_source_dir()
    #print source_dir
    f = PropertyFile()
    f.process(os.path.join(source_dir, "sidechain_angles.props"))
    for aa in f.properties['aminoacids'].split(","):
      #print aa + f.properties[aa+".chis"]
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

  def get_rotamers (self, residue_name) :
    rotamers = {}
    aa = residue_name.lower()
    if (not aa in self.rotamersForAA) :
      return None
    rotamer_list = self.rotamersForAA.get(aa)
    for rotamer in rotamer_list :
      if (rotamer != "") :
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
        if self.show_errors: print resName + " is an unknown residue type"
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

  def extract_chi_atoms (self, angleName, res, atom_dict=None) :
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
      if (testAtom != None) :
        angleAtoms.append(atomNamesMap.get(testAtom.name))
      else:
        return None
    return angleAtoms

def makeAtomDict(res):
  atomNamesMap = {}
  for atom in res.atoms():
    atomNamesMap[atom.name] = atom
  return atomNamesMap

def collect_sidechain_chi_angles (pdb_hierarchy, atom_selection=None) :
  angle_lookup = SidechainAngles(False)
  residue_chis = []
  if atom_selection is not None:
    if (isinstance(atom_selection, flex.bool)):
      actual_selection = atom_selection.iselection()
    elif (isinstance(atom_selection, flex.size_t)):
      actual_selection = atom_selection
  if atom_selection is None:
    actual_selection = flex.bool(
      len(pdb_hierarchy.atoms()),
      True).iselection()
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        for residue in conformer.residues() :
          n_chi = angle_lookup.chisPerAA.get(residue.resname.lower(), 0)
          try :
            n_chi = int(n_chi)
          except ValueError :
            continue
          chis = []
          altloc = residue.atoms()[0].fetch_labels().altloc
          i_seqs = []
          for i in range(1, n_chi+1) :
            atoms = angle_lookup.extract_chi_atoms("chi%d" % i, residue)
            if atoms is None:
              pass
            else :
              i_seqs = [ atom.i_seq for atom in atoms ]
              chis.append(group_args(chi_id=i, i_seqs=i_seqs))
          atoms_in_selection = True
          if atom_selection is not None:
            for i_seq in i_seqs:
              if i_seq not in actual_selection:
                atoms_in_selection = False
                break
          if (len(chis) > 0) and (atoms_in_selection) :
            residue_info = group_args(
              residue_name=residue.resname,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              chis=chis)
            residue_chis.append(residue_info)
  return residue_chis

def collect_residue_torsion_angles (pdb_hierarchy,
                                    atom_selection=None,
                                    chi_angles_only=False) :
  get_class = iotbx.pdb.common_residue_names_get_class
  residue_torsions = []

  ### chi angles ###
  residue_chis = collect_sidechain_chi_angles(
                   pdb_hierarchy=pdb_hierarchy,
                   atom_selection=atom_selection)
  residue_torsions = residue_chis
  ##################

  if atom_selection is not None:
    if (isinstance(atom_selection, flex.bool)):
      actual_selection = atom_selection.iselection()
    elif (isinstance(atom_selection, flex.size_t)):
      actual_selection = atom_selection
  if atom_selection is None:
    actual_selection = flex.bool(
      len(pdb_hierarchy.atoms()),
      True).iselection()
  previous_residue = None
  next_residue = None
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        for i_res, residue in enumerate(conformer.residues()) :
          if (get_class(residue.resname) != "common_amino_acid"):
            continue
          if i_res < (len(conformer.residues())-1):
            next_residue = conformer.residues()[i_res+1]
          else:
            next_residue = None
          torsions = []
          curN = None
          curCA = None
          curC = None
          for atom in residue.atoms():
            if atom.name == " N  ":
              curN = atom
            elif atom.name == " CA ":
              curCA = atom
            elif atom.name == " C  ":
              curC = atom
          if curN is None or curCA is None or curC is None:
            continue
          if not chi_angles_only:
            ### omega ###
            if previous_residue is not None:
              prevCA = None
              prevC = None
              for atom in previous_residue.atoms():
                if atom.name == " CA ":
                  prevCA = atom
                elif atom.name == " C  ":
                  prevC = atom
              if prevCA is not None and prevC is not None:
                atoms_in_selection = True
                if atom_selection is not None:
                  for atom in [prevCA, prevC, curN, curCA]:
                    if atom.i_seq not in actual_selection:
                      atoms_in_selection = False
                      break
                if atoms_in_selection:
                  omega = \
                    mmtbx.rotamer.omega_from_atoms(prevCA, prevC, curN, curCA)
                  if omega is not None:
                    i_seqs = [prevCA.i_seq,
                              prevC.i_seq,
                              curN.i_seq,
                              curCA.i_seq]
                    torsions.append(group_args(chi_id="omega", i_seqs=i_seqs))
            ###########

            ### phi ###
            if previous_residue is not None:
              prevC = None
              for atom in previous_residue.atoms():
                if atom.name == " C  ":
                  prevC = atom
              if prevC is not None:
                atoms_in_selection = True
                if atom_selection is not None:
                  for atom in [prevC, curN, curCA, curC]:
                    if atom.i_seq not in actual_selection:
                      atoms_in_selection = False
                      break
                if atoms_in_selection:
                  phi = mmtbx.rotamer.phi_from_atoms(prevC, curN, curCA, curC)
                  if phi is not None:
                    i_seqs = [prevC.i_seq,
                              curN.i_seq,
                              curCA.i_seq,
                              curC.i_seq]
                    torsions.append(group_args(chi_id="phi", i_seqs=i_seqs))
            ###########

            ### psi ###
            if next_residue is not None:
              nextN = None
              for atom in next_residue.atoms():
                if atom.name == " N  ":
                  nextN = atom
              if nextN is not None:
                atoms_in_selection = True
                if atom_selection is not None:
                  for atom in [curN, curCA, curC, nextN]:
                    if atom.i_seq not in actual_selection:
                      atoms_in_selection = False
                      break
                if atoms_in_selection:
                  psi = mmtbx.rotamer.psi_from_atoms(curN, curCA, curC, nextN)
                  if psi is not None:
                    i_seqs = [curN.i_seq,
                              curCA.i_seq,
                              curC.i_seq,
                              nextN.i_seq]
                    torsions.append(group_args(chi_id="psi", i_seqs=i_seqs))
            ###########

          ### c-beta ###
          #curCB = None
          #for atom in residue.atoms():
          #  if atom.name == " CB ":
          #    curCB = atom
          #if curCB is not None:
          #  atoms_in_selection = True
          #  for atom in [curN, curC, curCA, curCB]:
          #    if atom.i_seq not in atom_selection:
          #      atoms_in_selection = False
          #  if atoms_in_selection:
          #    ncab = \
          #      mmtbx.rotamer.improper_ncab_from_atoms(curN,
          #                                             curC,
          #                                             curCA,
          #                                             curCB)
          #    i_seqs = [curN.i_seq, curC.i_seq, curCA.i_seq, curCB.i_seq]
          #    torsions.append(group_args(chi_id="ncab", i_seqs=i_seqs))
          #    cnab = \
          #      mmtbx.rotamer.improper_ncab_from_atoms(curC,
          #                                             curN,
          #                                             curCA,
          #                                             curCB)
          #    i_seqs = [curC.i_seq,
          #              curN.i_seq,
          #              curCA.i_seq,
          #              curCB.i_seq]
          #    torsions.append(group_args(chi_id="cnab", i_seqs=i_seqs))
          ##############

          altloc = residue.atoms()[0].fetch_labels().altloc
          atoms_in_selection = True
          if (len(torsions) > 0) and (atoms_in_selection) :
            residue_info = group_args(
              residue_name=residue.resname,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              chis=torsions)
            residue_torsions.append(residue_info)
          previous_residue = residue
  return residue_torsions
