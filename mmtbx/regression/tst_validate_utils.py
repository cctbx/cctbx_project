from __future__ import division
#(jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.clashscore import clashscore
from mmtbx.validation.rna_validate import rna_validate
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from iotbx import pdb
from libtbx.test_utils import show_diff, Exception_expected
from libtbx.utils import Sorry
import libtbx.load_env

import sys, os

#{{{ exercise_rna_validate
def exercise_rna_validate():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb2goz_refmac_tls.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb2goz_refmac_tls.ent) not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)
  rv=rna_validate()
  rv.analyze_pdb(pdb_io=pdb_io)
  assert len(rv.pucker_outliers) == 2
  assert len(rv.bond_outliers) == 2
  assert len(rv.angle_outliers) == 0
  assert len(rv.suite_validation) == 3

#{{{ exercise_clashscore
def exercise_clashscore():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1zff.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb1zff.ent) not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)

  r = clashscore()
  score, bad_clashes = clashscore.analyze_clashes(r,pdb_io)
  clash_out = 'clashscore = %f' % score['']
  assert not show_diff(clash_out, "clashscore = 37.974684")
#}}}

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

  r = ramalyze()
  output, output_list = r.analyze_pdb(pdb_io=pdb_io, outliers_only=True)
  assert output.count("OUTLIER") == 100
  assert output.count("Favored") == 0
  assert output.count("Allowed") == 0
  assert output.count("General") == 64
  assert output.count("Glycine") == 6
  assert output.count("Trans-proline") == 1
  assert output.count("Cis-proline") == 0
  assert output.count("Pre-proline") == 4
  assert output.count("Isoleucine or valine") == 25

  output, output_list = r.analyze_pdb(pdb_io=pdb_io, outliers_only=False)
  assert output.count("OUTLIER") == 100
  assert output.count("Favored") == 461
  assert output.count("Allowed") == 162
  assert output.count("General") == 513
  assert output.count("Glycine") == 39
  assert output.count("Trans-proline") == 23
  assert output.count("Cis-proline") == 0
  assert output.count("Pre-proline") == 21
  assert output.count("Isoleucine or valine") == 127
  numtotal = r.get_phi_psi_residues_count()
  assert r.get_outliers_count_and_fraction()  == (100, 100./numtotal)
  assert r.get_allowed_count_and_fraction()   == (162, 162./numtotal)
  assert r.get_favored_count_and_fraction()   == (461, 461./numtotal)
  assert r.get_general_count_and_fraction()   == (513, 513./numtotal)
  assert r.get_gly_count_and_fraction()       == (39, 39./numtotal)
  assert r.get_trans_pro_count_and_fraction() == (23, 23./numtotal)
  assert r.get_cis_pro_count_and_fraction()   == (0, 0./numtotal)
  assert r.get_prepro_count_and_fraction()    == (21, 21./numtotal)
  assert r.get_ileval_count_and_fraction()    == (127, 127./numtotal)
  assert numtotal == 75+154+494
  output_lines = output.splitlines()
  assert len(output_lines) == 723
  assert output_lines[0] == "A  15  SER:35.07:-83.26:131.88:Favored:General"
  assert output_lines[1] == "A  16  SER:0.74:-111.53:71.36:Allowed:General"
  assert output_lines[168] == "A 191  ASP:2.66:-42.39:121.87:Favored:Pre-proline"
  assert output_lines[169] == "A 192  PRO:0.31:-39.12:-31.84:Allowed:Trans-proline"
  assert output_lines[713] == "B 368  LYS:56.44:-62.97:-53.28:Favored:General"
  assert output_lines[714] == "B 369  GLU:8.89:-44.36:-45.50:Favored:General"
  assert output_lines[715] == "B 370  LYS:40.00:-50.00:-39.06:Favored:General"
  assert output_lines[716] == "B 371  VAL:68.24:-60.38:-51.85:Favored:Isoleucine or valine"
  assert output_lines[717] == "B 372  LEU:0.02:-61.13:-170.23:OUTLIER:General"
  assert output_lines[718] == "B 373  ARG:0.02:60.09:-80.26:OUTLIER:General"
  assert output_lines[719] == "B 374  ALA:0.13:-37.21:-36.12:Allowed:General"
  assert output_lines[720] == "B 375  LEU:11.84:-89.81:-41.45:Favored:General"
  assert output_lines[721] == "B 376  ASN:84.33:-58.30:-41.39:Favored:General"
  assert output_lines[722] == "B 377  GLU:30.88:-56.79:-21.74:Favored:General"

  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_ramalyze(): input pdb (pdb1jxt.ent) not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)

  r = ramalyze()
  output, output_list = r.analyze_pdb(pdb_io=pdb_io, outliers_only=True)
  assert output.count("Favored") == 0
  assert output.count("Allowed") == 0
  assert output.count("OUTLIER") == 0

  output, output_list = r.analyze_pdb(pdb_io=pdb_io, outliers_only=False)
  assert output.count("Favored") == 47
  assert output.count("Allowed") == 1
  assert output.count("OUTLIER") == 0
  assert output.count("General") == 27
  assert output.count("Glycine") == 4
  assert output.count("Trans-proline") == 4
  assert output.count("Cis-proline") == 0
  assert output.count("Pre-proline") == 5
  assert output.count("Isoleucine or valine") == 8
  numtotal = r.get_phi_psi_residues_count()
  assert r.get_outliers_count_and_fraction()  == (0, 0./numtotal)
  assert r.get_allowed_count_and_fraction()   == (1, 1./numtotal)
  assert r.get_favored_count_and_fraction()   == (47, 47./numtotal)
  assert r.get_general_count_and_fraction()   == (27, 27./numtotal)
  assert r.get_gly_count_and_fraction()       == (4, 4./numtotal)
  assert r.get_trans_pro_count_and_fraction() == (4, 4./numtotal)
  assert r.get_cis_pro_count_and_fraction()   == (0, 0./numtotal)
  assert r.get_prepro_count_and_fraction()    == (5, 5./numtotal)
  assert r.get_ileval_count_and_fraction()    == (8, 8./numtotal)
  output_lines = output.splitlines()
  assert len(output_lines) == 48
  assert output_lines[0] == "A   2  ATHR:33.85:-106.92:144.23:Favored:General"
  assert output_lines[1] == "A   2  BTHR:37.07:-97.44:137.00:Favored:General"
  assert output_lines[6] == "A   7  AILE:98.76:-61.91:-44.35:Favored:Isoleucine or valine"
  assert output_lines[7] == "A   7  BILE:61.50:-56.21:-51.56:Favored:Isoleucine or valine"
  assert output_lines[8] == "A   8  AVAL:23.11:-50.35:-49.64:Favored:Isoleucine or valine"
  assert output_lines[9] == "A   8  BVAL:12.01:-83.20:-12.14:Favored:Isoleucine or valine"
  assert output_lines[10] == "A   8  CVAL:73.11:-61.22:-36.49:Favored:Isoleucine or valine"
  assert output_lines[44] == "A  43  AASP:51.81:-94.64:5.45:Favored:General"
  assert output_lines[45] == "A  43  BASP:56.98:-88.69:-0.12:Favored:General"
  assert output_lines[46] == "A  44  TYR:1.76:-133.10:58.75:Allowed:General"
  assert output_lines[47] == "A  45  ALA:57.37:-86.61:-8.57:Favored:General"
  # 2plx excerpt (unusual icode usage)
  pdb_io = pdb.input(source_info=None, lines="""\
ATOM   1468  N   GLY A 219       3.721  21.322  10.752  1.00 14.12           N
ATOM   1469  CA  GLY A 219       3.586  21.486  12.188  1.00 14.85           C
ATOM   1470  C   GLY A 219       4.462  20.538  12.995  1.00 15.63           C
ATOM   1471  O   GLY A 219       5.513  20.090  12.512  1.00 14.55           O
ATOM   1472  N   CYS A 220       4.036  20.213  14.235  1.00 15.02           N
ATOM   1473  CA  CYS A 220       4.776  19.228  15.068  1.00 15.56           C
ATOM   1474  C   CYS A 220       3.773  18.322  15.741  1.00 14.69           C
ATOM   1475  O   CYS A 220       2.799  18.828  16.338  1.00 15.54           O
ATOM   1476  CB  CYS A 220       5.620  19.906  16.174  1.00 15.72           C
ATOM   1477  SG  CYS A 220       6.762  21.133  15.448  1.00 15.45           S
ATOM   1478  N   ALA A 221A      4.054  17.017  15.707  1.00 14.77           N
ATOM   1479  CA  ALA A 221A      3.274  16.015  16.507  1.00 14.01           C
ATOM   1480  C   ALA A 221A      1.774  15.992  16.099  1.00 14.50           C
ATOM   1481  O   ALA A 221A      0.875  15.575  16.881  1.00 14.46           O
ATOM   1482  CB  ALA A 221A      3.440  16.318  17.935  1.00 12.28           C
ATOM   1483  N   GLN A 221       1.523  16.390  14.848  1.00 14.52           N
ATOM   1484  CA  GLN A 221       0.159  16.391  14.325  1.00 15.19           C
ATOM   1485  C   GLN A 221      -0.229  15.044  13.717  1.00 14.43           C
ATOM   1486  O   GLN A 221       0.641  14.280  13.307  1.00 16.88           O
ATOM   1487  CB  GLN A 221       0.002  17.491  13.272  1.00 16.41           C
ATOM   1488  CG  GLN A 221       0.253  18.906  13.805  1.00 16.52           C
ATOM   1489  CD  GLN A 221      -0.640  19.181  14.995  1.00 17.87           C
ATOM   1490  OE1 GLN A 221      -1.857  19.399  14.826  1.00 13.54           O
ATOM   1491  NE2 GLN A 221      -0.050  19.149  16.228  1.00 16.18           N
ATOM   1492  N   LYS A 222      -1.537  14.773  13.694  1.00 14.34           N
ATOM   1493  CA  LYS A 222      -2.053  13.536  13.125  1.00 15.07           C
ATOM   1494  C   LYS A 222      -1.679  13.455  11.655  1.00 14.88           C
ATOM   1495  O   LYS A 222      -1.856  14.424  10.883  1.00 14.32           O
""")
  output, output_list = r.analyze_pdb(pdb_io=pdb_io, outliers_only=False)
  assert (len(output_list) == 3)

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

  r = rotalyze()
  output, output_list = r.analyze_pdb(pdb_io, outliers_only=True)
  #print output
  assert output.count("OUTLIER") == 113
  #print output.count(":")
  assert output.count(":") == 791
  output_lines = output.splitlines()
  assert len(output_lines) == 113
  for lines in output_lines:
    assert float(lines[11:14]) <= 1.0

  output, output_list = r.analyze_pdb(pdb_io, outliers_only=False)
  assert output.count("OUTLIER") == 113
  assert output.count(":") == 4501
  assert output.count("p") == 121
  assert output.count("m") == 333
  assert output.count("t") == 495
  output_lines = output.splitlines()
  assert len(output_lines) == 643
  assert output_lines[0]   == "A  14  MET:1.00:3.3:29.2:173.3:287.9::ptm"
  assert output_lines[1]   == "A  15  SER:1.00:0.1:229.0::::OUTLIER"
  assert output_lines[2]   == "A  16  SER:1.00:4.2:277.9::::m"
  assert output_lines[42]  == "A  58  ASN:1.00:2.0:252.4:343.6:::m-20"
  assert output_lines[43]  == "A  59  ILE:1.00:2.0:84.2:186.7:::pt"
  assert output_lines[168] == "A 202  GLU:1.00:0.4:272.7:65.9:287.8::OUTLIER"
  assert output_lines[169] == "A 203  ILE:1.00:5.0:292.9:199.6:::mt"
  assert output_lines[450] == "B 154  THR:1.00:0.1:356.0::::OUTLIER"
  assert output_lines[587] == "B 316  TYR:1.00:5.4:153.7:68.6:::t80"
  assert output_lines[394] == "B  86  ASP:1.00:2.2:321.4:145.1:::m-20"
  assert output_lines[641] == "B 377  GLU:1.00:45.3:311.7:166.2:160.1::mt-10"
  assert output_lines[642] == "B 378  THR:1.00:23.5:309.4::::m"

  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_ramalyze(): input pdb (pdb1jxt.ent) not available"
    return
  pdb_io = pdb.input(file_name=regression_pdb)
  r = rotalyze()
  output, output_list = r.analyze_pdb(pdb_io, outliers_only=True)
  #print output
  assert output == "A  29  BTYR:0.35:0.3:191.3:322.7:::OUTLIER"

  output, output_list = r.analyze_pdb(pdb_io, outliers_only=False)
  assert not show_diff(output,"""\
A   1  THR:1.00:96.6:299.5::::m
A   2  ATHR:0.67:55.0:56.1::::p
A   2  BTHR:0.33:93.8:298.1::::m
A   3  CYS:1.00:28.5:310.5::::m
A   4  CYS:1.00:89.0:293.1::::m
A   5  PRO:1.00:90.6:30.2::::Cg_endo
A   6  SER:1.00:84.6:68.4::::p
A   7  AILE:0.45:62.7:290.8:178.2:::mt
A   7  BILE:0.55:14.3:284.4:298.4:::mm
A   8  AVAL:0.50:3.6:156.7::::t
A   8  BVAL:0.30:9.7:71.3::::p
A   8  CVAL:0.20:74.3:172.1::::t
A  10  AARG:0.65:23.4:176.8:66.5:63.9:180.0:tpp180
A  10  BARG:0.35:20.1:176.8:72.8:66.4:171.9:tpp180
A  11  SER:1.00:49.6:300.9::::m
A  12  AASN:0.50:96.0:286.1:343.8:::m-20
A  12  BASN:0.50:99.2:288.4:337.6:::m-20
A  13  APHE:0.65:42.3:187.2:276.4:::t80
A  13  BPHE:0.35:84.7:179.6:263.1:::t80
A  14  ASN:1.00:95.9:289.6:333.0:::m-20
A  15  VAL:1.00:47.4:168.2::::t
A  16  CYS:1.00:44.7:176.5::::t
A  17  ARG:1.00:23.6:289.7:282.8:288.6:158.7:mmm180
A  18  LEU:1.00:75.1:287.2:173.3:::mt
A  19  PRO:1.00:43.6:24.4::::Cg_endo
A  21  THR:1.00:8.5:314.0::::m
A  22  APRO:0.55:78.5:333.5::::Cg_exo
A  23  AGLU:0.50:92.5:290.9:187.1:341.8::mt-10
A  23  BGLU:0.50:94.5:292.0:183.8:339.2::mt-10
A  25  ALEU:0.50:96.7:294.4:173.6:::mt
A  26  CYS:1.00:92.2:295.0::::m
A  28  THR:1.00:37.5:52.9::::p
A  29  ATYR:0.65:23.0:161.8:67.8:::t80
A  29  BTYR:0.35:0.3:191.3:322.7:::OUTLIER
A  30  ATHR:0.70:68.5:57.4::::p
A  30  BTHR:0.30:8.8:78.1::::p
A  32  CYS:1.00:69.2:301.7::::m
A  33  ILE:1.00:37.5:66.5:173.4:::pt
A  34  AILE:0.70:66.6:303.6:167.6:::mt
A  34  BILE:0.30:33.9:308.5:296.8:::mm
A  35  ILE:1.00:48.4:62.4:170.0:::pt
A  36  PRO:1.00:36.1:22.5::::Cg_endo
A  39  ATHR:0.70:18.3:311.0::::m
A  39  BTHR:0.30:17.7:288.8::::m
A  40  CYS:1.00:99.0:294.4::::m
A  41  PRO:1.00:61.4:34.4::::Cg_endo
A  43  AASP:0.75:29.6:56.5:340.3:::p-10
A  43  BASP:0.25:45.3:59.6:349.3:::p-10
A  44  TYR:1.00:85.6:290.9:85.1:::m-85
A  46  ASN:1.00:34.0:301.6:117.9:::m120""")

  pdb_str = """\
ATOM   2527  N   LEU A 261     -31.022 -24.808 107.479  1.00 28.22           N
ATOM   2528  CA  LEU A 261     -30.054 -23.719 107.237  1.00 21.77           C
ATOM   2529  C   LEU A 261     -30.582 -22.773 106.168  1.00 27.64           C
ATOM   2530  O   LEU A 261     -29.841 -21.977 105.561  1.00 26.70           O
ATOM   2531  CB  LEU A 261     -28.696 -24.276 106.874  1.00 22.58           C
ATOM   2532  CG  LEU A 261     -28.135 -25.066 108.060  1.00 40.89           C
ATOM   2533  CD1 LEU A 261     -26.892 -25.858 107.664  1.00 46.72           C
ATOM   2534  CD2 LEU A 261     -27.806 -24.109 109.202  1.00 38.88           C
ATOM   2535  H   LEU A 261     -31.201 -25.277 106.781  1.00 33.87           H
ATOM   2536  HA  LEU A 261     -29.950 -23.204 108.064  1.00 26.12           H
ATOM   2537  HB2 LEU A 261     -28.781 -24.874 106.115  1.00 27.10           H
ATOM   2538  HB3 LEU A 261     -28.088 -23.548 106.670  1.00 27.10           H
ATOM   2539  HG  LEU A 261     -28.806 -25.693 108.373  1.00 49.07           H
ATOM   2540 HD11 LEU A 261     -26.570 -26.338 108.430  1.00 56.07           H
ATOM   2541 HD12 LEU A 261     -27.124 -26.473 106.965  1.00 56.07           H
ATOM   2542 HD13 LEU A 261     -26.219 -25.247 107.353  1.00 56.07           H
ATOM   2543 HD21 LEU A 261     -28.608 -23.653 109.468  1.00 46.66           H
ATOM   2544 HD22 LEU A 261     -27.455 -24.612 109.941  1.00 46.66           H
ATOM   2545 HD23 LEU A 261     -27.153 -23.474 108.899  1.00 46.66           H
ATOM   2546  N   GLY A 262     -31.887 -22.863 105.948  1.00 23.68           N
ATOM   2547  CA  GLY A 262     -32.572 -21.935 105.075  1.00 21.87      85   C
ATOM   2548  C   GLY A 262     -33.718 -22.620 104.386  1.00 27.32           C
ATOM   2549  O   GLY A 262     -33.943 -23.822 104.556  1.00 23.10           O
ATOM   2550  H   GLY A 262     -32.399 -23.459 106.298  1.00 28.42           H
ATOM   2551  HA2 GLY A 262     -32.916 -21.189 105.591  1.00 26.25      85   H
ATOM   2552  HA3 GLY A 262     -31.958 -21.598 104.405  1.00 26.25      85   H
ATOM   2553  N   SER A 263     -34.460 -21.830 103.628  1.00 24.62           N
ATOM   2554  CA  SER A 263     -35.631 -22.290 102.921  1.00 27.15           C
ATOM   2555  C   SER A 263     -35.594 -21.761 101.492  1.00 22.14           C
ATOM   2556  O   SER A 263     -34.723 -20.945 101.159  1.00 21.01           O
ATOM   2557  CB  SER A 263     -36.839 -21.713 103.619  1.00 25.73           C
ATOM   2558  OG  SER A 263     -36.907 -22.232 104.922  1.00 26.84           O
ATOM   2559  H   SER A 263     -34.296 -20.995 103.507  1.00 29.54           H
ATOM   2560  HA  SER A 263     -35.680 -23.269 102.917  1.00 32.58           H
ATOM   2561  HB2 SER A 263     -36.754 -20.747 103.661  1.00 30.87           H
ATOM   2562  HB3 SER A 263     -37.641 -21.960 103.132  1.00 30.87           H
ATOM   2563  HG  SER A 263     -37.560 -21.925 105.312  1.00 32.20           H
"""
  pdb_io = pdb.input(source_info=None, lines=pdb_str)
  try :
    rotalyze().analyze_pdb(pdb_io)
  except Sorry, e :
    assert ("GLY A 262" in str(e))
  else :
    raise Exception_expected
#}}}

def run():
  verbose = "--verbose" in sys.argv[1:]
  if (not libtbx.env.has_module(name="phenix")):
    print \
      "Skipping exercise_rna_validate():" \
      " phenix not available"
  else:
    exercise_rna_validate()
  exercise_ramalyze()
  exercise_rotalyze()
  if (not libtbx.env.has_module(name="probe")):
    print \
      "Skipping exercise_clashscore():" \
      " probe not available"
  else:
    exercise_clashscore()
  print "OK"

if (__name__ == "__main__"):
  run()
