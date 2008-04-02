# (jEdit options) :folding=explicit:collapseFolds=1:
import libtbx.load_env
import sys, os, getopt
try:
  from iotbx import pdb
  #from iotbx.pdb.hybrid_36 import hy36decode
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

#from mmtbx import monomer_library
#from mmtbx.monomer_library import server
#import mmtbx.monomer_library.pdb_interpretation
from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer import ramachandran_eval
from cctbx import geometry_restraints


#{{{ flag routines
#flag routines-----------------------------------------------------------------------------------
def help():
  version()
  print """USAGE:  mmtbx.ramalyze file.pdb

FLAGS:
  -h    Print this help message
  -v    Display version information
  -c    Print a change log
  -o    Only print outliers
"""
def changes():
  print "\nversion 0.1 080326 - First version\n"
def version():
  print "\nversion 0.1 080326 - Copyright 2007, Vincent Chen\n"
#------------------------------------------------------------------------------------------------
#}}}

#{{{ parse_cmdline
#parse the command line--------------------------------------------------------------------------
def parse_cmdline(params):
  try:
    opts, args = getopt.getopt( sys.argv[1:], 'hcvo',['help', 'changes', 'version', 'outliersonly'] )
  except getopt.GetoptError:
    help()
    sys.exit()
  for o, a in opts:
    if o in ("-h", "--help"):
      help()
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
    sys.exit(help())
  elif len(args) > 1:
    sys.stderr.write("\n**ERROR: too many input files specified\n")
    sys.exit(help())
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
  print analyze_pdb(pdb_io, params["outliersonly"])
#}}}

#{{{ analyze_pdb
def analyze_pdb(pdb_io, outliers_only=None):
  hierarchy      = pdb_io.construct_hierarchy()
  analysis = ""
  for model in hierarchy.models():
    for chain in model.chains():
      prevRes, prevC, resN, resCA, resO = None, None, None, None, None;
      for i_conformer, conformer in enumerate(chain.conformers()):
        #print str(conformer.residues()[-3].name)
        residues = conformer.residues()
        for i in range(len(residues)):
          residue = residues[i]
          phi = get_phi(residues, i)
          psi = get_psi(residues, i)

          resType = None
          if (residue.resname[0:3] == "GLY"):
            resType = "glycine"
          elif (residue.resname[0:3] == "PRO"):
            resType = "proline"
          elif (isPrePro(residues, i)):
            resType = "prepro"
          else:
            resType = "general"

          r = ramachandran_eval.RamachandranEval()
          value = 0
          if (phi is not None and psi is not None):
            value = r.evaluate(resType, [phi, psi])
            ramaType = evaluateScore(resType, value)
            #print params["outliersonly"]
            if (not outliers_only or isOutlier(resType, value)):
              analysis += '%s%4s %s:%.2f:%.2f:%.2f:%s:%s\n' % (chain.id,residue.resseq,residue.resname,value*100,phi,psi,ramaType,resType.capitalize())
              #print str(residue.seq).rjust(4,' ') + " " + residue.name + ":" + str(value) + ":" + str(phi) + ":" + str(psi) + ":OUTLIER:" + resType.capitalize()
              #print type(residue.seq)
  return analysis.rstrip()
#}}}

#{{{ get_phi
def get_phi(residues, i):
  prevC, resN, resCA, resC = None, None, None, None;
  if (i > 0 and i < len(residues)):
    prev = residues[i-1]
    if (prev is not None):
      for atom in prev.atoms():
        if (atom.name == " C  "): prevC = atom
    res = residues[i]
    if (res is not None):
      for atom in res.atoms():
        if (atom.name == " N  "): resN = atom
        if (atom.name == " CA "): resCA = atom
        if (atom.name == " C  "): resC = atom
    if (prevC is not None and resN is not None and resCA is not None and resC is not None):
      b = geometry_restraints.bond(
        sites=[prevC.xyz,resN.xyz],
        distance_ideal=1,
        weight=1)
      # check to see if residues are actually bonded.
      if (b.distance_model > 4): return None
      d = geometry_restraints.dihedral(
        sites=[prevC.xyz,resN.xyz,resCA.xyz,resC.xyz],
        angle_ideal=-40,
        weight=1)
      return d.angle_model
#}}}

#{{{ get_psi
def get_psi(residues, i):
  resN, resCA, resC, nextN = None, None, None, None;
  if (i >= 0 and i < len(residues) - 1):
    next = residues[i+1]
    if (next is not None):
      for atom in next.atoms():
        if (atom.name == " N  "): nextN = atom
    res = residues[i]
    if (res is not None):
      for atom in res.atoms():
        if (atom.name == " N  "): resN = atom
        if (atom.name == " CA "): resCA = atom
        if (atom.name == " C  "): resC = atom
    if (nextN is not None and resN is not None and resCA is not None and resC is not None):
      b = geometry_restraints.bond(
        sites=[resC.xyz,nextN.xyz],
        distance_ideal=1,
        weight=1)
      if (b.distance_model > 4): return None
      d = geometry_restraints.dihedral(
        sites=[resN.xyz,resCA.xyz,resC.xyz,nextN.xyz],
        angle_ideal=-40,
        weight=1)
      return d.angle_model
#}}}

#{{{ isPrePro
def isPrePro(residues, i):
  if (i < 0 or i >= len(residues) - 1): return False
  else:
    next = residues[i+1]
    if (next.resname[0:3] == "PRO"): return True
  return False
#}}}

#{{{ isOutlier
def isOutlier(resType, value):
  if (resType == "general"):
    if (value < 0.0005): return True
    else: return False
  else:
    if (value < 0.002): return True
    else: return False
#}}}

#{{{ evaluateScore
def evaluateScore(resType, value):
  if (value >= 0.02): return "Favored"
  if (resType == "general"):
    if (value >= 0.0005): return "Allowed"
    else: return "OUTLIER"
  else:
    if (value >= 0.0020): return "Allowed"
    else: return "OUTLIER"
#}}}

if __name__ == "__main__":
  params = {}
  params["outliersonly"]=False
  args = parse_cmdline(params)
  run(args, params)
