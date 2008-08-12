#(jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.command_line import ramalyze
from mmtbx.command_line import rotalyze
from mmtbx.command_line import cbetadev
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
import libtbx.load_env

import sys, os, getopt

#{{{ exercise_cbetadev
def exercise_cbetadev():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb1jxt.ent) not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)

  output = cbetadev.analyze_pdb(regression_pdb,pdb_io,True)
  assert not show_diff(output, """\
pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:
pdb1jxt :a:ile: A:   7 :  0.260: -46.47:   0.45:a:
pdb1jxt :b:val: A:   8 :  0.258:  80.92:   0.30:b:
pdb1jxt :c:val: A:   8 :  0.641: -53.98:   0.20:c:
pdb1jxt :b:thr: A:  30 :  0.812: -76.98:   0.30:b:
pdb1jxt :b:thr: A:  39 :  0.924:  56.41:   0.30:b:
pdb1jxt :b:asp: A:  43 :  0.500:   7.56:   0.25:b:""")

  output = cbetadev.analyze_pdb(regression_pdb,pdb_io, False)
  assert not show_diff(output, """\
pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:
pdb1jxt : :thr: A:   1 :  0.102:  11.27:   1.00: :
pdb1jxt :a:thr: A:   2 :  0.022: -49.31:   0.67:a:
pdb1jxt : :cys: A:   3 :  0.038: 103.68:   1.00: :
pdb1jxt : :cys: A:   4 :  0.047:-120.73:   1.00: :
pdb1jxt : :pro: A:   5 :  0.069:-121.41:   1.00: :
pdb1jxt : :ser: A:   6 :  0.052: 112.87:   1.00: :
pdb1jxt :a:ile: A:   7 :  0.260: -46.47:   0.45:a:
pdb1jxt :b:ile: A:   7 :  0.153: 122.97:   0.55:b:
pdb1jxt :a:val: A:   8 :  0.184:-155.36:   0.50:a:
pdb1jxt :b:val: A:   8 :  0.258:  80.92:   0.30:b:
pdb1jxt :c:val: A:   8 :  0.641: -53.98:   0.20:c:
pdb1jxt : :ala: A:   9 :  0.061: -82.84:   1.00: :
pdb1jxt :a:arg: A:  10 :  0.023: 172.24:   1.00:a:
pdb1jxt : :ser: A:  11 :  0.028:-129.11:   1.00: :
pdb1jxt :a:asn: A:  12 :  0.021: -80.80:   0.50:a:
pdb1jxt :b:asn: A:  12 :  0.199:  50.01:   0.50:b:
pdb1jxt :a:phe: A:  13 :  0.067: -37.32:   0.65:a:
pdb1jxt :b:phe: A:  13 :  0.138:  19.24:   0.35:b:
pdb1jxt : :asn: A:  14 :  0.065: -96.35:   1.00: :
pdb1jxt : :val: A:  15 :  0.138: -96.63:   1.00: :
pdb1jxt : :cys: A:  16 :  0.102: -28.64:   1.00: :
pdb1jxt : :arg: A:  17 :  0.053:-106.79:   1.00: :
pdb1jxt : :leu: A:  18 :  0.053:-141.51:   1.00: :
pdb1jxt : :pro: A:  19 :  0.065:-146.95:   1.00: :
pdb1jxt : :thr: A:  21 :  0.086:  53.80:   1.00: :
pdb1jxt :a:pro: A:  22 :  0.092: -83.39:   0.55:a:
pdb1jxt :a:glu: A:  23 :  0.014:-179.53:   0.50:a:
pdb1jxt :b:glu: A:  23 :  0.050:-179.78:   0.50:b:
pdb1jxt : :ala: A:  24 :  0.056: -88.96:   1.00: :
pdb1jxt : :leu: A:  25 :  0.084:-106.42:   1.00: :
pdb1jxt : :cys: A:  26 :  0.074: -94.70:   1.00: :
pdb1jxt : :ala: A:  27 :  0.056: -62.15:   1.00: :
pdb1jxt : :thr: A:  28 :  0.056:-114.82:   1.00: :
pdb1jxt :a:tyr: A:  29 :  0.068:   0.22:   0.65:a:
pdb1jxt :a:thr: A:  30 :  0.180: 103.27:   0.70:a:
pdb1jxt :b:thr: A:  30 :  0.812: -76.98:   0.30:b:
pdb1jxt : :cys: A:  32 :  0.029: -84.07:   1.00: :
pdb1jxt : :ile: A:  33 :  0.048:-119.17:   1.00: :
pdb1jxt : :ile: A:  34 :  0.045:  99.02:   1.00: :
pdb1jxt : :ile: A:  35 :  0.052:-128.24:   1.00: :
pdb1jxt : :pro: A:  36 :  0.084:-142.29:   1.00: :
pdb1jxt : :ala: A:  38 :  0.039:  50.02:   1.00: :
pdb1jxt :a:thr: A:  39 :  0.093: -96.63:   0.70:a:
pdb1jxt :b:thr: A:  39 :  0.924:  56.41:   0.30:b:
pdb1jxt : :cys: A:  40 :  0.013:-144.12:   1.00: :
pdb1jxt : :pro: A:  41 :  0.039: -97.09:   1.00: :
pdb1jxt :a:asp: A:  43 :  0.130:-146.91:   0.75:a:
pdb1jxt :b:asp: A:  43 :  0.500:   7.56:   0.25:b:
pdb1jxt : :tyr: A:  44 :  0.085:-143.63:   1.00: :
pdb1jxt : :ala: A:  45 :  0.055:  33.32:   1.00: :
pdb1jxt : :asn: A:  46 :  0.066: -50.46:   1.00: :""")

#{{{ exercise_ramalyze
def exercise_ramalyze():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_ramalyze(): input pdb (jcm.pdb) not available"
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print "Skipping exercise_ramalyze(): rotarama_data directory not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)

  output = ramalyze.analyze_pdb(pdb_io, True)
  assert output.count("OUTLIER") == 75
  assert output.count("Favored") == 0
  assert output.count("Allowed") == 0
  assert output.count("General") == 65
  assert output.count("Glycine") == 5
  assert output.count("Proline") == 1
  assert output.count("Prepro") == 4

  output = ramalyze.analyze_pdb(pdb_io, False)
  assert output.count("OUTLIER") == 75
  assert output.count("Favored") == 496
  assert output.count("Allowed") == 154
  assert output.count("General") == 642
  assert output.count("Glycine") == 39
  assert output.count("Proline") == 23
  assert output.count("Prepro") == 21
  output_lines = output.splitlines()
  assert len(output_lines) == 725
  assert output_lines[0] == "A  15 SER:39.85:-83.26:131.88:Favored:General"
  assert output_lines[1] == "A  16 SER:0.93:-111.53:71.36:Allowed:General"
  assert output_lines[168] == "A 191 ASP:2.90:-42.39:121.87:Favored:Prepro"
  assert output_lines[169] == "A 192 PRO:3.65:-39.12:-31.84:Favored:Proline"
  assert output_lines[715] == "B 368 LYS:62.62:-62.97:-53.28:Favored:General"
  assert output_lines[716] == "B 369 GLU:9.58:-44.36:-45.50:Favored:General"
  assert output_lines[717] == "B 370 LYS:37.37:-50.00:-39.06:Favored:General"
  assert output_lines[718] == "B 371 VAL:71.48:-60.38:-51.85:Favored:General"
  assert output_lines[719] == "B 372 LEU:0.04:-61.13:-170.23:OUTLIER:General"
  assert output_lines[720] == "B 373 ARG:0.03:60.09:-80.26:OUTLIER:General"
  assert output_lines[721] == "B 374 ALA:0.57:-37.21:-36.12:Allowed:General"
  assert output_lines[722] == "B 375 LEU:13.45:-89.81:-41.45:Favored:General"
  assert output_lines[723] == "B 376 ASN:84.52:-58.30:-41.39:Favored:General"
  assert output_lines[724] == "B 377 GLU:32.22:-56.79:-21.74:Favored:General"
#}}}

#{{{ exercise_rotalyze
def exercise_rotalyze():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_rotalyze(): input pdb (jcm.pdb) not available"
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print "Skipping exercise_rotalyze(): rotarama_data directory not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)

  output = rotalyze.analyze_pdb(pdb_io, True)
  assert output.count("OUTLIER") == 113
  assert output.count(":") == 678
  output_lines = output.splitlines()
  assert len(output_lines) == 113
  for lines in output_lines:
    assert float(lines[10:13]) <= 1.0

  output = rotalyze.analyze_pdb(pdb_io, False)
  assert output.count("OUTLIER") == 113
  assert output.count(":") == 3858
  assert output.count("p") == 120
  assert output.count("m") == 333
  assert output.count("t") == 496
  output_lines = output.splitlines()
  assert len(output_lines) == 643
  assert output_lines[0]   == "A  14 MET:3.3:29.2:173.3:287.9::ptm"
  assert output_lines[1]   == "A  15 SER:0.1:229.0::::OUTLIER"
  assert output_lines[2]   == "A  16 SER:4.2:277.9::::m"
  assert output_lines[42]  == "A  58 ASN:2.0:252.4:343.6:::m-20"
  assert output_lines[43]  == "A  59 ILE:2.0:84.2:186.7:::pt"
  assert output_lines[168] == "A 202 GLU:0.4:272.7:65.9:107.8::OUTLIER"
  assert output_lines[169] == "A 203 ILE:5.0:292.9:199.6:::mt"
  assert output_lines[450] == "B 154 THR:0.1:356.0::::OUTLIER"
  assert output_lines[587] == "B 316 TYR:5.4:153.7:68.6:::t80"
  assert output_lines[394] == "B  86 ASP:2.2:321.4:145.1:::m-20"
  assert output_lines[641] == "B 377 GLU:45.3:311.7:166.2:160.1::mt-10"
  assert output_lines[642] == "B 378 THR:23.5:309.4::::m"
#}}}

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_ramalyze()
  exercise_rotalyze()
  exercise_cbetadev()
  print "OK"

if (__name__ == "__main__"):
  run()
