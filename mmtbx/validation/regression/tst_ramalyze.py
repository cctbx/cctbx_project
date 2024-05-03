
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import show_diff, approx_equal
import libtbx.load_env
from libtbx.easy_pickle import loads, dumps
from six.moves import cStringIO as StringIO
import os.path
from mmtbx.validation import ramalyze
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long

import time
import json

def exercise_ramalyze():
  from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
  import iotbx.pdb
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_ramalyze(): input pdb (jcm.pdb) not available")
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print("Skipping exercise_ramalyze(): rotarama_data directory not available")
    return
  # Exercise 1
  pdb_in = iotbx.pdb.input(file_name=regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  pdb_io = iotbx.pdb.input(file_name=regression_pdb)
  hierarchy.atoms().reset_i_seq()
  r = ramalyze.ramalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  out = StringIO()
  r.show_old_output(out=out)
  output = out.getvalue()
  assert output.count("OUTLIER") == 100
  assert output.count("Favored") == 0
  assert output.count("Allowed") == 0
  assert output.count("General") == 64
  assert output.count("Glycine") == 6
  assert output.count("Trans-proline") == 1
  assert output.count("Cis-proline") == 0
  assert output.count("Pre-proline") == 4
  assert output.count("Isoleucine or valine") == 25
  assert (len(r.outlier_selection()) == 494)
  outlier_ids = set([])
  atoms = hierarchy.atoms()
  for i_seq in r.outlier_selection():
    atom = atoms[i_seq]
    atom_group = atoms[i_seq].parent()
    outlier_ids.add(atom_group.id_str())
  outliers1 = sorted([ o.atom_group_id_str() for o in r.results ])
  outliers2 = sorted(list(outlier_ids))
  assert (outliers1 == outliers2)

  r = ramalyze.ramalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  for unpickle in [False, True] :
    if unpickle :
      r = loads(dumps(r))
    for outlier in r.results :
      assert (len(outlier.xyz) == 3)
    out = StringIO()
    r.show_old_output(out=out, verbose=False)
    output = out.getvalue()
    assert output.count("OUTLIER") == 100
    assert output.count("Favored") == 463
    assert output.count("Allowed") == 162
    assert output.count("General") == 514
    assert output.count("Glycine") == 39
    assert output.count("Trans-proline") == 23
    assert output.count("Cis-proline") == 0
    assert output.count("Pre-proline") == 21
    assert output.count("Isoleucine or valine") == 128
    numtotal = r.get_phi_psi_residues_count()
    assert r.get_outliers_count_and_fraction()  == (100, 100./numtotal)
    assert r.get_allowed_count_and_fraction()   == (162, 162./numtotal)
    assert r.get_favored_count_and_fraction()   == (463, 463./numtotal)
    assert r.get_general_count_and_fraction()   == (514, 514./numtotal)
    assert r.get_gly_count_and_fraction()       == (39, 39./numtotal)
    assert r.get_trans_pro_count_and_fraction() == (23, 23./numtotal)
    assert r.get_cis_pro_count_and_fraction()   == (0, 0./numtotal)
    assert r.get_prepro_count_and_fraction()    == (21, 21./numtotal)
    assert r.get_ileval_count_and_fraction()    == (128, 128./numtotal)
    #assert numtotal == 75+154+494 #reasons for this math unclear
    assert numtotal == 725
    output_lines = output.splitlines()
    assert len(output_lines) == 725
    selected_lines = []
    for x in [0, 1, 168, 169, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724]:
      selected_lines.append(output_lines[x])
    assert not show_diff("\n".join(selected_lines), """\
 A  15  SER:35.07:-83.26:131.88:Favored:General
 A  16  SER:0.74:-111.53:71.36:Allowed:General
 A 191  ASP:2.66:-42.39:121.87:Favored:Pre-proline
 A 192  PRO:0.31:-39.12:-31.84:Allowed:Trans-proline
 B 368  LYS:56.44:-62.97:-53.28:Favored:General
 B 369  GLU:8.89:-44.36:-45.50:Favored:General
 B 370  LYS:40.00:-50.00:-39.06:Favored:General
 B 371  VAL:68.24:-60.38:-51.85:Favored:Isoleucine or valine
 B 372  LEU:0.02:-61.13:-170.23:OUTLIER:General
 B 373  ARG:0.02:60.09:-80.26:OUTLIER:General
 B 374  ALA:0.13:-37.21:-36.12:Allowed:General
 B 375  LEU:11.84:-89.81:-41.45:Favored:General
 B 376  ASN:84.33:-58.30:-41.39:Favored:General
 B 377  GLU:30.88:-56.79:-21.74:Favored:General""")
    assert (len(r.outlier_selection()) == 494)

  # Exercise 2
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  pdb_in = iotbx.pdb.input(file_name=regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  r = ramalyze.ramalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  out = StringIO()
  r.show_old_output(out=out)
  output = out.getvalue()
  assert output.count("Favored") == 0
  assert output.count("Allowed") == 0
  assert output.count("OUTLIER") == 0
  r = ramalyze.ramalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  for unpickle in [False, True] :
    if unpickle :
      r = loads(dumps(r))
    out = StringIO()
    r.show_old_output(out=out, verbose=False)
    output = out.getvalue()
    assert output.count("Favored") == 50
    assert output.count("Allowed") == 1
    assert output.count("OUTLIER") == 0
    assert output.count("General") == 29
    assert output.count("Glycine") == 4
    assert output.count("Trans-proline") == 5
    assert output.count("Cis-proline") == 0
    assert output.count("Pre-proline") == 5
    assert output.count("Isoleucine or valine") == 8
    numtotal = r.get_phi_psi_residues_count()
    assert r.get_outliers_count_and_fraction()  == (0, 0./numtotal)
    assert r.get_allowed_count_and_fraction()   == (1, 1./numtotal)
    assert r.get_favored_count_and_fraction()   == (43, 43./numtotal)
    #print r.get_general_count_and_fraction()
    assert r.get_general_count_and_fraction()   == (25, 25./numtotal)
    assert r.get_gly_count_and_fraction()       == (4, 4./numtotal)
    assert r.get_trans_pro_count_and_fraction() == (5, 5./numtotal)
    assert r.get_cis_pro_count_and_fraction()   == (0, 0./numtotal)
    assert r.get_prepro_count_and_fraction()    == (5, 5./numtotal)
    assert r.get_ileval_count_and_fraction()    == (5, 5./numtotal)
    output_lines = output.splitlines()
    assert len(output_lines) == 51
    selected_lines = []
    for x in [0, 1, 5, 6, 7, 8, 9, 47, 48, 49, 50]:
      selected_lines.append(output_lines[x])
    assert not show_diff("\n".join(selected_lines), """\
 A   2 ATHR:33.85:-106.92:144.23:Favored:General
 A   3 ACYS:47.07:-132.54:137.26:Favored:General
 A   7 AILE:98.76:-61.91:-44.35:Favored:Isoleucine or valine
 A   7 BILE:61.50:-56.21:-51.56:Favored:Isoleucine or valine
 A   8 AVAL:23.11:-50.35:-49.64:Favored:Isoleucine or valine
 A   8 BVAL:12.01:-83.20:-12.14:Favored:Isoleucine or valine
 A   8 CVAL:73.11:-61.22:-36.49:Favored:Isoleucine or valine
 A  43 AASP:51.81:-94.64:5.45:Favored:General
 A  43 BASP:56.98:-88.69:-0.12:Favored:General
 A  44  TYR:1.76:-133.10:58.75:Allowed:General
 A  45  ALA:57.37:-86.61:-8.57:Favored:General""")

  # Exercise 3: 2plx excerpt (unusual icode usage)
  import iotbx.pdb
  hierarchy = iotbx.pdb.input(source_info=None, lines="""\
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
""").construct_hierarchy()
  r = ramalyze.ramalyze(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  assert (len(r.results) == 3)

def exercise_favored_regions():
  assert ramalyze.get_favored_regions(0) == [(-99, 119), (-63, -43), (53, 43), (60,-120)]

def exercise_constants():
  #
  # if this test fails, somebody changed ramalyze constants. The same constants
  # are declared also in mmtbx/validation/ramachandran/rama8000_tables.h
  # It is essential to keep both places consistent.
  #
  assert ramalyze.res_types == ["general", "glycine", "cis-proline",
      "trans-proline", "pre-proline", "isoleucine or valine"]
  assert ramalyze.RAMA_GENERAL == 0
  assert ramalyze.RAMA_GLYCINE == 1
  assert ramalyze.RAMA_CISPRO == 2
  assert ramalyze.RAMA_TRANSPRO == 3
  assert ramalyze.RAMA_PREPRO == 4
  assert ramalyze.RAMA_ILE_VAL == 5
  assert ramalyze.RAMALYZE_OUTLIER == 0
  assert ramalyze.RAMALYZE_ALLOWED == 1
  assert ramalyze.RAMALYZE_FAVORED == 2
  assert ramalyze.RAMALYZE_ANY == 3
  assert ramalyze.RAMALYZE_NOT_FAVORED == 4

def exercise_ramalyze_json(test_mmcif=False):
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/jcm.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_ramalyze(): input pdb (jcm.pdb) not available")
    return
  if (find_rotarama_data_dir(optional=True) is None):
    print("Skipping exercise_ramalyze(): rotarama_data directory not available")
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
  ramalyze_json = ramalyze.ramalyze(pdb_hierarchy=m.get_hierarchy(), outliers_only=True).as_JSON()
  rmjson_dict = json.loads(ramalyze_json)
  #import pprint
  #pprint.pprint(rmjson_dict)
  assert len(rmjson_dict['flat_results'])==100, "tst_ramalyze json output not returning correct number of values"
  assert approx_equal(rmjson_dict['flat_results'][0]['phi'], 50.51521639791719), "tst_ramalyze json output first calculated phi dihedral angle not matching previous value"
  assert approx_equal(rmjson_dict['flat_results'][0]['psi'], -80.04604513007598), "tst_ramalyze json output first calculated psi dihedral angle not matching previous value"
  assert rmjson_dict['flat_results'][0]['rama_type']=='OUTLIER', "tst_ramalyze json output first rama_type not matching previous value"
  assert approx_equal(rmjson_dict['flat_results'][99]['phi'], 60.09378543010022), "tst_ramalyze json output last calculated phi dihedral angle not matching previous value"
  assert approx_equal(rmjson_dict['flat_results'][99]['psi'], -80.26327714086905), "tst_ramalyze json output last calculated psi dihedral angle not matching previous value"
  assert rmjson_dict['flat_results'][99]['rama_type']=='OUTLIER', "tst_ramalyze json output last rama_type not matching previous value"
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(rmjson_dict['hierarchical_results'], "OUTLIER")==100, "tst_ramalyze json hierarchical output total number of rama outliers changed"
  assert rmjson_dict['summary_results'][""]['num_allowed'] == 162, "tst_ramalyze json output summary total num_allowed not matching previous value"
  assert rmjson_dict['summary_results'][""]['num_favored'] == 463, "tst_ramalyze json output summary total num_favored not matching previous value"
  assert rmjson_dict['summary_results'][""]['num_outliers'] == 100, "tst_ramalyze json output summary total num_outliers not matching previous value"
  assert rmjson_dict['summary_results'][""]['num_residues'] == 725, "tst_ramalyze json output summary total num_residues not matching previous value"
  return rmjson_dict

if (__name__ == "__main__"):
  t0=time.time()
  exercise_ramalyze()
  exercise_favored_regions()
  exercise_constants()
  rm_dict = exercise_ramalyze_json()
  rm_dict_cif = exercise_ramalyze_json(test_mmcif=True)
  assert rm_dict['summary_results'] == rm_dict_cif['summary_results'], "tst_ramalyze summary results changed between pdb and cif version"
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
