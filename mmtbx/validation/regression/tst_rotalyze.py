
from __future__ import absolute_import, division, print_function
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from mmtbx.validation import rotalyze
from iotbx import pdb
from libtbx.test_utils import show_diff, Exception_expected, approx_equal
from libtbx.utils import Sorry
import libtbx.load_env
from libtbx.easy_pickle import loads, dumps
from six.moves import cStringIO as StringIO
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long
import os.path
import json
from six.moves import zip

def exercise_rotalyze():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_rotalyze(): input pdb (jcm.pdb) not available")
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print("Skipping exercise_rotalyze(): rotarama_data directory not available")
    return
  pdb_in = pdb.input(file_name=regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  pdb_io = pdb.input(file_name=regression_pdb)
  r = rotalyze.rotalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  out = StringIO()
  r.show_old_output(out=out, verbose=False)
  output = out.getvalue()
  assert output.count("OUTLIER") == 246, output.count("OUTLIER")
  assert output.count(":") == 984, output.count(":")
  output_lines = output.splitlines()
  assert len(output_lines) == 123
  for lines in output_lines:
    assert float(lines[12:15]) <= 1.0

  r = rotalyze.rotalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  for unpickle in [False, True] :
    if unpickle :
      r = loads(dumps(r))
    out = StringIO()
    r.show_old_output(out=out, verbose=False)
    for outlier in r.results :
      assert (len(outlier.xyz) == 3)
    output = out.getvalue()
    assert output.count("OUTLIER") == 246
    assert output.count(":") == 5144, output.count(":")
    assert output.count("p") == 120
    assert output.count("m") == 324
    assert output.count("t") == 486
    output_lines = output.splitlines()
    #for line in output_lines:
    #  print line
    #STOP()
    assert len(output_lines) == 643
    line_indices = [0,1,2,42,43,168,169,450,587,394,641,642]

#    top500 version
    line_values = [
     " A  14  MET:1.00:3.3:29.2:173.3:287.9::Favored:ptm",
     " A  15  SER:1.00:0.1:229.0::::OUTLIER:OUTLIER",
     " A  16  SER:1.00:4.2:277.9::::Favored:m",
     " A  58  ASN:1.00:2.0:252.4:343.6:::Favored:m-20",
     " A  59  ILE:1.00:2.0:84.2:186.7:::Allowed:pt",
     " A 202  GLU:1.00:0.4:272.7:65.9:287.8::OUTLIER:OUTLIER",
     " A 203  ILE:1.00:5.0:292.9:199.6:::Favored:mt",
     " B 154  THR:1.00:0.1:356.0::::OUTLIER:OUTLIER",
     " B 316  TYR:1.00:5.4:153.7:68.6:::Favored:t80",
     " B  86  ASP:1.00:2.2:321.4:145.1:::Favored:m-20",
     " B 377  GLU:1.00:45.3:311.7:166.2:160.1::Favored:mt-10",
     " B 378  THR:1.00:23.5:309.4::::Favored:m"]
#    top8000 version
    line_values = [
     " A  14  MET:1.00:1.3:29.2:173.3:287.9::Allowed:ptm",
     " A  15  SER:1.00:0.1:229.0::::OUTLIER:OUTLIER",
     " A  16  SER:1.00:3.0:277.9::::Favored:m",
     " A  58  ASN:1.00:1.0:252.4:343.6:::Allowed:m-40",
     " A  59  ILE:1.00:0.5:84.2:186.7:::Allowed:pt",
     " A 202  GLU:1.00:0.0:272.7:65.9:287.8::OUTLIER:OUTLIER",
     " A 203  ILE:1.00:1.0:292.9:199.6:::Allowed:mt",
     " B 154  THR:1.00:0.0:356.0::::OUTLIER:OUTLIER",
     " B 316  TYR:1.00:4.1:153.7:68.6:::Favored:t80",
     " B  86  ASP:1.00:0.4:321.4:145.1:::Allowed:m-30",
     " B 377  GLU:1.00:15.0:311.7:166.2:160.1::Favored:mt-10",
     " B 378  THR:1.00:17.0:309.4::::Favored:m",
    ]
    for idx, val in zip(line_indices, line_values):
      assert (output_lines[idx] == val), (idx, output_lines[idx])

  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_ramalyze(): input pdb (pdb1jxt.ent) not available")
    return
  pdb_in = pdb.input(file_name=regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  pdb_io = pdb.input(file_name=regression_pdb)
  r = rotalyze.rotalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  out = StringIO()
  r.show_old_output(out=out, verbose=False)
  output = out.getvalue().strip()
  assert output == ""

  r = rotalyze.rotalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  for unpickle in [False, True] :
    if unpickle :
      r = loads(dumps(r))
    out = StringIO()
    r.show_old_output(out=out, verbose=False)
    output = out.getvalue()
    assert not show_diff(output,"""\
 A   1  THR:1.00:95.4:299.5::::Favored:m
 A   2 ATHR:0.67:49.5:56.1::::Favored:p
 A   2 BTHR:0.33:90.4:298.1::::Favored:m
 A   3  CYS:1.00:12.9:310.5::::Favored:m
 A   4  CYS:1.00:91.6:293.1::::Favored:m
 A   5  PRO:1.00:78.8:30.2:319.7:33.8::Favored:Cg_endo
 A   6  SER:1.00:90.1:68.4::::Favored:p
 A   7 AILE:0.45:49.6:290.8:178.2:::Favored:mt
 A   7 BILE:0.55:6.5:284.4:298.4:::Favored:mm
 A   8 AVAL:0.50:1.1:156.7::::Allowed:t
 A   8 BVAL:0.30:5.1:71.3::::Favored:p
 A   8 CVAL:0.20:69.8:172.1::::Favored:t
 A  10 AARG:0.65:24.7:176.8:66.5:63.9:180.0:Favored:tpp-160
 A  10 BARG:0.35:17.5:176.8:72.8:66.4:171.9:Favored:tpp-160
 A  11  SER:1.00:51.6:300.9::::Favored:m
 A  12 AASN:0.50:93.9:286.1:343.8:::Favored:m-40
 A  12 BASN:0.50:98.9:288.4:337.6:::Favored:m-40
 A  13 APHE:0.65:45.1:187.2:276.4:::Favored:t80
 A  13 BPHE:0.35:86.1:179.6:263.1:::Favored:t80
 A  14  ASN:1.00:95.2:289.6:333.0:::Favored:m-40
 A  15  VAL:1.00:42.3:168.2::::Favored:t
 A  16  CYS:1.00:40.8:176.5::::Favored:t
 A  17  ARG:1.00:21.4:289.7:282.8:288.6:158.7:Favored:mmm160
 A  18  LEU:1.00:65.0:287.2:173.3:::Favored:mt
 A  19  PRO:1.00:43.6:24.4:324.8:31.6::Favored:Cg_endo
 A  21  THR:1.00:5.7:314.0::::Favored:m
 A  22 APRO:0.55:87.5:333.5:34.0:333.8::Favored:Cg_exo
 A  23 AGLU:0.50:86.9:290.9:187.1:341.8::Favored:mt-10
 A  23 BGLU:0.50:91.7:292.0:183.8:339.2::Favored:mt-10
 A  25 ALEU:0.50:95.7:294.4:173.6:::Favored:mt
 A  26  CYS:1.00:83.0:295.0::::Favored:m
 A  28  THR:1.00:29.6:52.9::::Favored:p
 A  29 ATYR:0.65:18.5:161.8:67.8:::Favored:t80
 A  29 BTYR:0.35:0.4:191.3:322.7:::Allowed:t80
 A  30 ATHR:0.70:60.8:57.4::::Favored:p
 A  30 BTHR:0.30:6.6:78.1::::Favored:p
 A  32  CYS:1.00:61.4:301.7::::Favored:m
 A  33  ILE:1.00:36.6:66.5:173.4:::Favored:pt
 A  34 AILE:0.70:60.9:303.6:167.6:::Favored:mt
 A  34 BILE:0.30:31.4:308.5:296.8:::Favored:mm
 A  35  ILE:1.00:45.6:62.4:170.0:::Favored:pt
 A  36  PRO:1.00:36.2:22.5:330.5:24.8::Favored:Cg_endo
 A  39 ATHR:0.70:14.0:311.0::::Favored:m
 A  39 BTHR:0.30:13.1:288.8::::Favored:m
 A  40  CYS:1.00:81.4:294.4::::Favored:m
 A  41  PRO:1.00:35.4:34.4:317.5:33.1::Favored:Cg_endo
 A  43 AASP:0.75:24.8:56.5:340.3:::Favored:p0
 A  43 BASP:0.25:43.2:59.6:349.3:::Favored:p0
 A  44  TYR:1.00:85.3:290.9:85.1:::Favored:m-80
 A  46  ASN:1.00:38.7:301.6:117.9:::Favored:m110
""")

def exercise_2():
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

  pdb_str2 = """
ATOM    453  N   PRO A  47       8.633   6.370   5.022  1.00 13.79           N
ATOM    454  CA  PRO A  47       7.915   7.571   5.496  1.00 14.61           C
ATOM    455  C   PRO A  47       7.612   7.481   6.994  1.00 15.06           C
ATOM    456  O   PRO A  47       7.289   6.377   7.439  1.00 14.39           O
ATOM    457  CB  PRO A  47       6.639   7.559   4.651  1.00 16.24           C
ATOM    458  CG  PRO A  47       7.089   6.901   3.338  1.00 15.52           C
ATOM    459  CD  PRO A  47       7.990   5.773   3.833  1.00 14.40           C
ATOM    460  N   MSE A  48       7.754   8.528   7.779  1.00 15.13           N
ATOM    461  CA  MSE A  48       7.482   8.456   9.201  1.00 16.17           C
ATOM    462  C   MSE A  48       6.040   8.750   9.517  1.00 15.23           C
ATOM    463  O   MSE A  48       5.417   9.418   8.735  1.00 14.77           O
ATOM    464  CB  MSE A  48       8.165   9.538  10.023  1.00 19.62           C
ATOM    465  CG  MSE A  48       9.630   9.466  10.238  1.00 21.70           C
ATOM    466 SE   MSE A  48      10.022  10.161  12.050  0.70 37.95          SE
ATOM    467  CE  MSE A  48      11.268   8.720  12.235  1.00 28.72           C
ATOM    468  N   LYS A  49       5.519   8.291  10.645  1.00 13.93           N
ATOM    469  CA  LYS A  49       4.167   8.624  11.045  1.00 13.79           C
ATOM    470  C   LYS A  49       4.022  10.138  11.202  1.00 14.66           C
ATOM    471  O   LYS A  49       5.011  10.853  11.351  1.00 15.69           O
ATOM    472  CB  LYS A  49       3.797   7.915  12.349  1.00 13.33           C
ATOM    473  CG  LYS A  49       3.593   6.416  12.204  1.00 14.35           C
ATOM    474  CD  LYS A  49       2.121   6.071  12.044  1.00 16.45           C
ATOM    475  CE  LYS A  49       1.571   5.402  13.292  1.00 18.19           C
ATOM    476  NZ  LYS A  49       0.899   4.110  12.980  1.00 19.97           N
"""

  pdb_io = pdb.input(source_info=None, lines=pdb_str)
  hierarchy = pdb_io.construct_hierarchy()
  try :
    rotalyze.rotalyze(pdb_hierarchy=hierarchy)
  except Sorry as e :
    assert ("GLY A 262" in str(e))
  else :
    raise Exception_expected

  pdb_io = pdb.input(source_info=None, lines=pdb_str2)
  hierarchy = pdb_io.construct_hierarchy()
  r = rotalyze.rotalyze(pdb_hierarchy=hierarchy)
  out = StringIO()
  r.show_old_output(out=out, verbose=False)
  output = out.getvalue()
  assert output == """\
 A  47  PRO:1.00:86.4:329.3:41.3:324.9::Favored:Cg_exo
 A  48  MSE:0.70:0.3:287.6:214.8:138.3::OUTLIER:OUTLIER
 A  49  LYS:1.00:0.1:288.6:263.2:251.7:233.0:OUTLIER:OUTLIER
""", output

  r = rotalyze.rotalyze(pdb_hierarchy=hierarchy,
    data_version="8000")
  out = StringIO()
  r.show_old_output(out=out, verbose=False)
  assert (out.getvalue() == """\
 A  47  PRO:1.00:86.4:329.3:41.3:324.9::Favored:Cg_exo
 A  48  MSE:0.70:0.3:287.6:214.8:138.3::OUTLIER:OUTLIER
 A  49  LYS:1.00:0.1:288.6:263.2:251.7:233.0:OUTLIER:OUTLIER
"""), out.getvalue()

  try :
    r = rotalyze.rotalyze(pdb_hierarchy=hierarchy,
      data_version="9000")
  except ValueError :
    pass
  else :
    raise Exception_expected

  from mmtbx.rotamer.rotamer_eval import RotamerEval
  rotamer_manager = RotamerEval()
  results = []
  for model in hierarchy.models():
    for chain in model.chains():
      for residue in chain.residues():
        cur_rot = rotamer_manager.evaluate_residue(residue)
        results.append(cur_rot)
  assert results == ['Cg_exo', 'OUTLIER', 'OUTLIER']

def exercise_rotalyze_json(test_mmcif=False):
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_rotalyze(): input pdb (jcm.pdb) not available")
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print("Skipping exercise_rotalyze(): rotarama_data directory not available")
    return
  dm = DataManager()
  if test_mmcif:
    with open(regression_pdb) as f:
      pdb_jcm_str = f.read()
    pdb_jcm_str = convert_string_to_cif_long(pdb_jcm_str, chain_addition="LONGCHAIN")
    dm.process_model_str("1", pdb_jcm_str)
    m = dm.get_model("1")
  else:
    m = dm.get_model(regression_pdb)
  rotalyze_json = rotalyze.rotalyze(pdb_hierarchy=m.get_hierarchy(), outliers_only=True).as_JSON()
  rtjson_dict = json.loads(rotalyze_json)
  #import pprint
  #pprint.pprint(rtjson_dict)
  assert len(rtjson_dict['flat_results'])==123, "tst_rotalyze json output not returning correct number of values"
  assert approx_equal(rtjson_dict['flat_results'][0]['chi_angles'][0], 229.02299329063914), "tst_rotalyze json output first calculated chi dihedral angle not matching previous value"
  assert rtjson_dict['flat_results'][0]['rotamer_name']=='OUTLIER', "tst_rotalyze json output first rotamer_name not matching previous value"
  assert approx_equal(rtjson_dict['flat_results'][122]['chi_angles'][0], 328.0085051658891), "tst_rotalyze json output last calculated first chi dihedral angle not matching previous value"
  assert approx_equal(rtjson_dict['flat_results'][122]['chi_angles'][1], 352.23811983072466), "tst_rotalyze json output last calculated second chi dihedral angle not matching previous value"
  assert rtjson_dict['flat_results'][122]['rotamer_name']=='OUTLIER', "tst_rotalyze json output last rotamer_name not matching previous value"
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(rtjson_dict['hierarchical_results'], "OUTLIER")==246, "tst_rotalyze json hierarchical output total number of rota outliers changed"
  assert rtjson_dict['summary_results'][""]['num_allowed'] == 116, "tst_rotalyze json output summary total num_allowed not matching previous value"
  assert rtjson_dict['summary_results'][""]['num_favored'] == 404, "tst_rotalyze json output summary total num_favored not matching previous value"
  assert rtjson_dict['summary_results'][""]['num_outliers'] == 123, "tst_rotalyze json output summary total num_outliers not matching previous value"
  assert rtjson_dict['summary_results'][""]['num_residues'] == 643, "tst_rotalyze json output summary total num_residues not matching previous value"
  return rtjson_dict

if (__name__ == "__main__"):
  exercise_rotalyze()
  exercise_2()
  rt_dict = exercise_rotalyze_json()
  rt_dict_cif = exercise_rotalyze_json(test_mmcif=True)
  assert rt_dict['summary_results'] == rt_dict_cif['summary_results'], "tst_rotalyze summary results changed between pdb and cif version"

  print("OK")
