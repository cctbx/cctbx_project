import libtbx.load_env
from cctbx import geometry_restraints
import sys, os

def find_source_dir(optional=False):
  #result = libtbx.env.find_in_repositories("mmtbx")
  #if result is None:
  result = libtbx.env.find_in_repositories(
    os.path.join("mmtbx", "mmtbx", "rotamer"))
  if result is None and not optional:
    raise Sorry("""\
Can't seem to find mmtbx/mmtbx/rotamer/ directory.
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
  #rotamersForAA = {}

  def __init__(self):
    source_dir = find_source_dir()
    #print source_dir
    f = PropertyFile()
    f.process(os.path.join(source_dir, "sidechain_angles.props"))
    for aa in f.properties['aminoacids'].split(","):
      #print aa + f.properties[aa+".chis"]
      self.chisPerAA[aa] = f.properties[aa+".chis"] #gives aaName -> # of chis
      #print f.properties[aa+".angles"].split(",")
      anglelist = f.properties[aa+".angles"].split(",")
      #print anglelist.count('')
      self.anglesForAA[aa] = anglelist #aaName -> [mobile angles]
      for angle in anglelist:
        if angle != '':
          key = aa+"."+angle
          self.atomsForAngle[key] = f.properties[key].split(",") #aaName.angle -> atoms

  def measureChiAngles(self, res):
    resName = res.resname.lower().strip()
    try:
      numChis = int(self.chisPerAA[resName])
      values = []
      for i in range(numChis):
        values.append(self.measureAngle("chi"+str(i+1), res))
      return values
    except KeyError:
      resName + " is unknown"


  def measureAngle(self, angleName, res):
    atomNamesMap = makeAtomDict(res)
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
      if testAtom != None: angleAtoms.append(atomNamesMap.get(testAtom.name))
      else: raise AttributeError, "some sidechain atoms are missing!"
      #testAtom = res.atoms[res.atoms().index(angleAtom)]
    d = geometry_restraints.dihedral(
      sites=[angleAtoms[0].xyz,angleAtoms[1].xyz,angleAtoms[2].xyz,angleAtoms[3].xyz],
      angle_ideal=-40,
      weight=1)
    return d.angle_model

def makeAtomDict(res):
  atomNamesMap = {}
  for atom in res.atoms():
    atomNamesMap[atom.name] = atom
  return atomNamesMap
