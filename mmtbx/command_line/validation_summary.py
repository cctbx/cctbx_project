
from libtbx.utils import Sorry, Usage
import cStringIO
import os
import sys

def run (args, out=sys.stdout) :
  if (len(args) == 0) :
    raise Usage("""
mmtbx.validation_summary model.pdb

Prints a brief summary of validation criteria, including Ramachandran
statistics, rotamer outliers, clashscore, C-beta deviations, plus R-factors
and RMS(bonds)/RMS(angles) if found in PDB header.  (This is primarily used
for evaluating the output of refinement tests; general users are advised to
run phenix.model_vs_data or the validation GUI.)
""")
  pdb_file = args[0]
  if (not os.path.isfile(pdb_file)) :
    raise Sorry("Not a file: %s" % pdb_file)
  from mmtbx.validation import ramalyze, rotalyze, cbetadev, clashscore
  log = cStringIO.StringIO()
  cs = clashscore.clashscore()
  clash_score = cs.run(args=[pdb_file], out=log, quiet=True)['']
  rama = ramalyze.ramalyze()
  rama.run(args=[pdb_file], out=log, quiet=True)
  rama_fav = rama.fav_percent
  rama_out = rama.out_percent
  rota = rotalyze.rotalyze()
  rota.run(args=[pdb_file], out=log, quiet=True)
  rota_out = rota.out_percent
  cbeta = cbetadev.cbetadev()
  cbeta_out = len(cbeta.run(args=[pdb_file, "cbetadev.outliers_only=True"],
    out=log, quiet=True))
  pdb_lines = open(pdb_file, "r").readlines()
  r_work = None
  r_free = None
  rms_bonds = None
  rms_angles = None
  for line in pdb_lines :
    if (line.startswith("REMARK   3")) :
      if ("Final:" in line) :
        fields = line.split()
        for i, field in enumerate(fields) :
          if (field == "r_work") :
            r_work = float(fields[i+2])
          elif (field == "r_free") :
            r_free = float(fields[i+2])
          elif (field == "bonds") :
            rms_bonds = float(fields[i+2])
          elif (field == "angles") :
            rms_angles = float(fields[i+2])
        break
      elif ("R VALUE            (WORKING SET)" in line) :
        r_work = float(line.split(":")[1].strip())
      elif ("FREE R VALUE                    " in line) :
        r_free = float(line.split(":")[1].strip())
    elif (line.startswith("REMARK 200")) :
      break
  print >> out, ""
  print >> out, "Validation summary for %s:" % pdb_file
  print >> out, "  Ramachandran outliers = %6.2f %%" % rama_out
  print >> out, "               favored  = %6.2f %%" % rama_fav
  print >> out, "  Rotamer outliers      = %6.2f %%" % rota_out
  print >> out, "  C-beta deviations     = %6d" % cbeta_out
  print >> out, "  Clashscore            = %6.2f" % clash_score
  if (r_work is not None) :
    print >> out, "  R-work                = %8.4f" % r_work
  if (r_free is not None) :
    print >> out, "  R-free                = %8.4f" % r_free
  if (rms_bonds is not None) :
    print >> out, "  RMS(bonds)            = %8.4f" % rms_bonds
  if (rms_angles is not None) :
    print >> out, "  RMS(angles)           = %6.2f" % rms_angles
  print >> out, ""

if (__name__ == "__main__") :
  run(sys.argv[1:])
