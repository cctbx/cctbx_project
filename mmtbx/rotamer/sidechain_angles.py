import scitbx.math
from libtbx.utils import Sorry
import libtbx.load_env
from libtbx import group_args
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

  def __init__(self, show_errs):
    self.show_errors = show_errs
    source_dir = find_source_dir()
    #print source_dir
    f = PropertyFile()
    f.process(os.path.join(source_dir, "sidechain_angles.props"))
    for aa in f.properties['aminoacids'].split(","):
      #print aa + f.properties[aa+".chis"]
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

  def measureChiAngles(self, res, atom_dict = None):
    resName = res.resname.lower().strip()
    try:
      numChis = int(self.chisPerAA[resName])
      values = []
      for i in range(numChis):
        values.append(self.measureAngle("chi"+str(i+1), res, atom_dict))
      return values
    except KeyError:
      if self.show_errors: print resName + " is an unknown residue type"

#  def measureChiAnglesDict(self, atom_dict, resName):
#    try:
#      numChis = int(self.chisPerAA[resName])
#      values = []
#      for i in range(numChis):
#        values.append(self.measureAngle("chi"+str(i+1), res))
#      return values
#    except KeyError:
#      resName + " is unknown"

  def measureAngle(self, *args, **kwds):
    angleAtoms = self.extract_chi_atoms(*args, **kwds)
    return scitbx.math.dihedral_angle(
      sites=[a.xyz for a in angleAtoms], deg=True)

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
      #print testAtom.name + res.name
      if (testAtom != None) :
        angleAtoms.append(atomNamesMap.get(testAtom.name))
      else:
        raise AttributeError("Missing sidechain atoms!")
    return angleAtoms

def makeAtomDict(res):
  atomNamesMap = {}
  for atom in res.atoms():
    atomNamesMap[atom.name] = atom
  return atomNamesMap

def collect_sidechain_chi_angles (pdb_hierarchy, atom_selection=None) :
  angle_lookup = SidechainAngles(False)
  residue_chis = []
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
          for i in range(1, n_chi+1) :
            try :
              atoms = angle_lookup.extract_chi_atoms("chi%d" % i, residue)
            except RuntimeError, e :
              pass
            else :
              i_seqs = [ atom.i_seq for atom in atoms ]
              chis.append(group_args(chi_id=i, i_seqs=i_seqs))
          if (len(chis) > 0) :
            residue_info = group_args(
              residue_name=residue.resname,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              chis=chis)
            residue_chis.append(residue_info)
  return residue_chis
