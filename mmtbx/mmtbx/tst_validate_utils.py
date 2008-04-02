#(jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.command_line import ramalyze
from mmtbx.command_line import rotalyze
from iotbx import pdb
#import iotbx.pdb.interpretation
from cctbx.array_family import flex
import libtbx.load_env

import sys, os, getopt

def find_rotarama_data_dir(optional=False):
  result = libtbx.env.find_in_repositories("rotarama_data")
  if result is None:
    result = libtbx.env.find_in_repositories(
      os.path.join("ext_ref_files", "rotarama_data"))
    if result is None and not optional:
      raise Sorry("""\
Can't find ext_ref_files/rotarama_data/ directory:
  Please run
    svn co svn://quiddity.biochem.duke.edu:21/phenix/rotarama_data
  to resolve this problem.""")
  return result

#{{{ exercise_ramalyze
def exercise_ramalyze():
  #pdb_io = iotbx.pdb.input(source_info=None, lines=flex.split_lines("""\
#ATOM   1748  CD  PRO A 260      75.791 136.003  84.241  1.00 37.59           C
#"""))
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
  print "OK"

if (__name__ == "__main__"):
  run()
