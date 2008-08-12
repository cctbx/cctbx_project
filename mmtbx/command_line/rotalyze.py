#(jEdit options) :folding=explicit:collapseFolds=1:
import libtbx.load_env
import sys, os, getopt
#try:
from iotbx import pdb
  #from iotbx.pdb.hybrid_36 import hy36decode
#except ImportError, e:
#  print "iotbx not loaded"
#  sys.exit()

from cctbx import geometry_restraints
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer.rotamer_eval import RotamerID


#{{{ flag routines
#flag routines-----------------------------------------------------------------------------------
def help():
  version()
  print """USAGE:  mmtbx.rotalyze file.pdb

FLAGS:
  -h    Print this help message
  -v    Display version information
  -c    Print a change log
  -o    Only print outliers
  -e    Show errors
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
    opts, args = getopt.getopt( sys.argv[1:], 'hcvoe',['help', 'changes', 'version', 'outliersonly', 'errors'] )
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
    if o in ("-e", "--errors"):
      params["showerrors"] = True

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
  #mon_lib_srv = server.server()
  filename = args[0]
  #ener_lib = server.ener_lib()
  if os.path.exists(filename):
    pdb_io = pdb.input(filename)
  elif filename.find("\n")>-1:
    pdb_io = pdb.input(source_info=None,
                       lines=flex.split_lines(filename))
  else:
    assert 0
  print analyze_pdb(pdb_io, params["outliersonly"], params["showerrors"])
#}}}

#{{{ analyze_pdb
def analyze_pdb(pdb_io, outliers_only=False, show_errors = False):
  sa = SidechainAngles()
  hierarchy      = pdb_io.construct_hierarchy()
  analysis = ""
  rot_id = rotamer_eval.RotamerID() # loads in the rotamer names
  #print rot_id.names
  for model in hierarchy.models():
    for chain in model.chains():
      for i_conformer, conformer in enumerate(chain.conformers()):
        #print str(conformer.residues()[-3].name)
        residues = conformer.residues()
        for i in range(len(residues)):
          residue = residues[i]
          try:
            chis = sa.measureChiAngles(residue)
            r = rotamer_eval.RotamerEval()
            if (chis is not None):
              value = r.evaluate(residue.resname.lower().strip(), chis)
              #print chain.id + str(residue.seq) + residue.name.strip() + str(value)
              if value != None:
                s = '%s%4s %s:%.1f' % (chain.id,residue.resseq,residue.resname.strip(),value*100)
                #wrap_chis = []
                #for i in range(len(chis)):
                #  wrap_chis.append(chis[i] % 360)
                #  if wrap_chis[i] < 0: wrap_chis[i] += 360
                wrap_chis = rot_id.wrap_chis(residue.resname.strip(), chis)
                for i in range(4):
                  s += ':'
                  if i < len(wrap_chis): s += '%.1f' % (wrap_chis[i])
                s += ':'
                if value < 0.01:
                  s += "OUTLIER\n"
                  if outliers_only: analysis += s
                else:
                  s += rot_id.identify(residue.resname, wrap_chis) + "\n"
                if not outliers_only: analysis += s
          except AttributeError:
            if show_errors: print '%s%4s %s is missing some sidechain atoms' % (chain.id, residue.resseq, residue.resname)
  return analysis.rstrip()
#}}}

if __name__ == "__main__":
  params = {}
  params["outliersonly"]=False
  params["showerrors"]=False
  args = parse_cmdline(params)
  run(args, params)
