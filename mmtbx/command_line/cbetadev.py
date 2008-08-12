import libtbx.load_env
import sys, os, getopt, math
try:
  from iotbx import pdb
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer import ramachandran_eval
from cctbx import geometry_restraints

#{{{ local_help
#flag routines-----------------------------------------------------------------------------------
def local_help():
  version()
  print """USAGE:  mmtbx.cbetadev file.pdb

FLAGS:
  -h    Print this help message
  -v    Display version information
  -c    Print a change log
  -o    Only print outliers (deviation >= 0.25)
"""
def changes():
  print "\nversion 0.10 080520 - First version\n"
def version():
  print "\nversion 0.10 080520 - Copyright 2008, Jeffrey J. Headd\n"
#------------------------------------------------------------------------------------------------
#}}}

#{{{ parse_cmdline
#parse the command line--------------------------------------------------------------------------
def parse_cmdline(params):
  try:
    opts, args = getopt.getopt( sys.argv[1:], 'hcvo',['help', 'changes', 'version', 'outliersonly'] )
  except getopt.GetoptError:
    local_help()
    sys.exit()
  for o, a in opts:
    if o in ("-h", "--help"):
      local_help()
      sys.exit()
    if o in ("-c", "--changes"):
      changes()
      sys.exit()
    if o in ("-v", "--version"):
      version()
      sys.exit()
    if o in ("-o", "--outliersonly"):
      params["outliersonly"] = True
  if len(args) < 1:
    sys.stderr.write("\n**ERROR: User must specify input PDB file\n")
    sys.exit(local_help())
  elif len(args) > 1:
    sys.stderr.write("\n**ERROR: too many input files specified\n")
    sys.exit(local_help())
  else:
    return args
#------------------------------------------------------------------------------------------------
#}}}

#{{{ run
def run(args, params=None, log=None):
  if (log is None): log = sys.stdout
  filename = args[0]
  if os.path.exists(filename):
    pdb_io = pdb.input(filename)
  elif filename.find("\n")>-1:
    pdb_io = pdb.input(source_info=None,
                       lines=flex.split_lines(filename))
  else:
    assert 0
  print analyze_pdb(filename, pdb_io, params["outliersonly"])
#}}}

#{{{ analyze_pdb
def analyze_pdb(filename, pdb_io, outliers_only=None):
  relevant_atom_names = {
    " CA ": None, " N  ": None, " C  ": None, " CB ": None} # FUTURE: set
  analysis = 'pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:\n'
  hierarchy = pdb_io.construct_hierarchy()
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for i_cf,cf in enumerate(rg.conformers()):
          for i_residue,residue in enumerate(cf.residues()):
            is_first = (i_cf == 0)
            is_alt_conf = False
            relevant_atoms = {}
            for atom in residue.atoms():
              if (atom.name in relevant_atom_names):
                relevant_atoms[atom.name] = atom
                if (len(atom.parent().altloc) != 0):
                  is_alt_conf = True
            if ((is_first or is_alt_conf) and len(relevant_atoms) == 4):
              resCA = relevant_atoms[" CA "]
              resN  = relevant_atoms[" N  "]
              resC  = relevant_atoms[" C  "]
              resCB = relevant_atoms[" CB "]
              dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal=idealized_calpha_angles(residue=residue)
              betaNCAB = construct_fourth(resN,resCA,resC,dist,angleCAB,dihedralNCAB,method="NCAB")
              betaCNAB = construct_fourth(resN,resCA,resC,dist,angleNAB,dihedralCNAB,method="CNAB")
              betaxyz = [(betaNCAB[0]+betaCNAB[0])/2,(betaNCAB[1]+betaCNAB[1])/2,(betaNCAB[2]+betaCNAB[2])/2]
              betadist = distance(resCA.xyz,betaxyz)
              if(betadist != dist):
                distTemp = [(betaxyz[0]-resCA.xyz[0]),(betaxyz[1]-resCA.xyz[1]),(betaxyz[2]-resCA.xyz[2])]
                betaxyz = [(resCA.xyz[0]+distTemp[0]*dist/betadist),(resCA.xyz[1]+distTemp[1]*dist/betadist),(resCA.xyz[2]+distTemp[2]*dist/betadist)]
              if(residue.resname != "GLY"):
                dev = distance(resCB.xyz,betaxyz)
                if(dev >=0.25 or outliers_only==False):
                  d = geometry_restraints.dihedral(sites=[resN.xyz,resCA.xyz,betaxyz,resCB.xyz],angle_ideal=0,weight=1)
                  dihedralNABB = d.angle_model

                  #internal version of dihedral calculation used to test difference
                  #dihedralTemp = dihedral4pt(resN.xyz,resCA.xyz,betaxyz,resCB.xyz)
                  #dihedralNABB = dihedralTemp

                  PDBfileStr = os.path.basename(filename)[:-4]
                  if (is_alt_conf):
                    altchar = cf.altloc.lower()
                  else:
                    altchar = " "
                  res=residue.resname.lower()
                  sub=chain.id
                  if(len(sub)==1):
                    sub=" "+sub
                  resnum=residue.resid()
                  resins=" "
                  occ = resCB.occ
                  analysis += '%s :%s:%s:%s:%4d%c:%7.3f:%7.2f:%7.2f:%s:\n' % (PDBfileStr,altchar,res,sub,int(resnum),resins,dev,dihedralNABB,occ,altchar)
  return analysis.rstrip()
#}}}

#{{{ idealized_calpha_angles
def idealized_calpha_angles(residue):
  if(residue.resname == "ALA"):
    dist = 1.536
    angleCAB = 110.1
    dihedralNCAB = 122.9
    angleNAB = 110.6
    dihedralCNAB = -122.6
    angleideal = 111.2
  elif(residue.resname == "PRO"):
    dist = 1.530
    angleCAB = 112.2
    dihedralNCAB = 115.1
    angleNAB = 103.0
    dihedralCNAB = -120.7
    angleideal = 111.8
  elif(residue.resname == "VAL") or (residue.resname == "THR") or (residue.resname == "ILE"):
    dist = 1.540
    angleCAB = 109.1
    dihedralNCAB = 123.4
    angleNAB = 111.5
    dihedralCNAB = -122.0
    angleideal = 111.2
  elif(residue.resname == "GLY"):
    dist = 1.10
    angleCAB = 109.3
    dihedralNCAB = 121.6
    angleNAB = 109.3
    dihedralCNAB = -121.6
    angleideal = 112.5
  else:
    dist = 1.530
    angleCAB = 110.1
    dihedralNCAB = 122.8
    angleNAB = 110.5
    dihedralCNAB = -122.6
    angleideal = 111.2
  return dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal
#}}}

#{{{ construct_fourth
def construct_fourth(resN,resCA,resC,dist,angle,dihedral,method="NCAB"):
  if (resN is not None and resCA is not None and resC is not None):
    if (method is "NCAB"):
      res0 = resN
      res1 = resC
      res2 = resCA
    elif (method is "CNAB"):
      res0 = resC
      res1 = resN
      res2 = resCA
    a = [(res2.xyz[0] - res1.xyz[0]),(res2.xyz[1] - res1.xyz[1]),(res2.xyz[2] - res1.xyz[2])]
    b = [(res0.xyz[0] - res1.xyz[0]),(res0.xyz[1] - res1.xyz[1]),(res0.xyz[2] - res1.xyz[2])]
    c = cross(a,b)
    cmag = mag(c)
    if(cmag > 0.000001):
      c=[(c[0]*dist/cmag),(c[1]*dist/cmag),(c[2]*dist/cmag)]
    c=[(c[0] + res2.xyz[0]),(c[1] + res2.xyz[1]),(c[2] + res2.xyz[2])]
    d = c
    angledhdrl = dihedral - 90
    a = [res1.xyz[0],res1.xyz[1],res1.xyz[2]]
    b = [res2.xyz[0],res2.xyz[1],res2.xyz[2]]
    newD = doaxisrot(d,angledhdrl,a,b)
    a = [(newD[0]-res2.xyz[0]),(newD[1]-res2.xyz[1]),(newD[2]-res2.xyz[2])]
    b = [(res1.xyz[0]-res2.xyz[0]),(res1.xyz[1]-res2.xyz[1]),(res1.xyz[2]-res2.xyz[2])]
    c = cross(a,b)
    cmag = mag(c)
    if(cmag > 0.000001):
      c=[(c[0]*dist/cmag),(c[1]*dist/cmag),(c[2]*dist/cmag)]
    angledhdrl = 90 - angle;
    a = [res2.xyz[0],res2.xyz[1],res2.xyz[2]]
    c = [(c[0]+res2.xyz[0]),(c[1]+res2.xyz[1]),(c[2]+res2.xyz[2])]
    b = c
    return doaxisrot(newD,angledhdrl,a,b)
#}}}

#{{{ doaxisrot
def doaxisrot(d,theta,res1,res2):
  LOK = True
  xx = res2[0] - res1[0]
  yy = res2[1] - res1[1]
  zz = res2[2] - res1[2]
  if(xx**2 + yy**2 + zz**2 > 0.000001):
    cosn1 = xx/(xx**2 + yy**2 + zz**2)**.5
    cosn2 = yy/(xx**2 + yy**2 + zz**2)**.5
    cosn3 = zz/(xx**2 + yy**2 + zz**2)**.5
  else:
    cosn1 = xx
    cosn2 = yy
    cosn3 = zz
  while(theta > 360):
    theta = theta - 360
  while(theta < -360):
    theta = theta + 360
  if(theta > 180):
    theta = theta - 360
  if(theta < -180):
    theta = theta + 360
  if(theta > -.0001 and theta < 0.0001):
    LOK = False
  if(LOK):
    costheta = math.cos(theta*3.14159*2/360)
    sintheta = math.sin(theta*3.14159*2/360)
    a1 = [(cosn1**2 + (1-cosn1**2)*costheta),(cosn1*cosn2*(1-costheta)+cosn3*sintheta),(cosn1*cosn3*(1-costheta) - cosn2*sintheta)]
    a2 = [(cosn1*cosn2*(1-costheta) - cosn3*sintheta),(cosn2**2 + (1-cosn2**2)*costheta),(cosn2*cosn3*(1-costheta) + cosn1*sintheta)]
    a3 = [(cosn1*cosn3*(1-costheta)+cosn2*sintheta),(cosn2*cosn3*(1-costheta)-cosn1*sintheta),(cosn3**2+(1-cosn3**2)*costheta)]
    f1 = [(d[0]-res2[0]),(d[1]-res2[1]),(d[2]-res2[2])]
    f2 = [(f1[0]*a1[0] + f1[1]*a2[0] + f1[2]*a3[0]),(f1[0]*a1[1] + f1[1]*a2[1] + f1[2]*a3[1]),(f1[0]*a1[2]+f1[1]*a2[2]+f1[2]*a3[2])]
    return [(f2[0]+res2[0]),(f2[1]+res2[1]),(f2[2]+res2[2])]
#}}}

#{{{ dihedral4pt
def dihedral4pt(p1,p2,p3,p4):
  a = [(p1[0]-p2[0]),(p1[1]-p2[1]),(p1[2]-p2[2])]
  b = [(p3[0]-p2[0]),(p3[1]-p2[1]),(p3[2]-p2[2])]
  d = cross(a,b)
  b = [(p2[0]-p3[0]),(p2[1]-p3[1]),(p2[2]-p3[2])]
  c = [(p4[0]-p3[0]),(p4[1]-p3[1]),(p4[2]-p3[2])]
  e = cross(b,c)
  dot = dotProduct(d,e)
  dmag = mag(d)
  emag = mag(e)
  if(dmag*emag < 0.0001):
    angle = 0.0
  else:
    angle = math.acos( dot/(dmag*emag) )
  anglehdrl = angle*360.0/(2*math.pi)
  f = cross(d,b)
  dot = dotProduct(f,e)
  fmag = mag(f)
  if(fmag*emag < 0.0001):
    angle = 0.0
  else:
    angle = math.acos( dot/(fmag*emag) )
  angle = angle*360.0/(2*math.pi)
  if(angle > 90.0):
    anglehdrl = -anglehdrl
  return anglehdrl
#}}}

#{{{ mag
def mag(v):
  return (v[0]**2 + v[1]**2 + v[2]**2)**0.5
#}}}

#{{{ cross
def cross(v1,v2):
  newcoords = [0,0,0]
  newcoords[0] = v1[1]*v2[2]-v1[2]*v2[1]
  newcoords[1] = v1[2]*v2[0]-v1[0]*v2[2]
  newcoords[2] = v1[0]*v2[1]-v1[1]*v2[0]
  return newcoords
#}}}

#{{{ dotProduct
def dotProduct(v1,v2):
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#}}}

#{{{ distance
def distance(a,b):
  return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
#}}}

if __name__ == "__main__":
  params = {}
  params["outliersonly"]=False
  args = parse_cmdline(params)
  run(args, params)
