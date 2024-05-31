
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
from libtbx.easy_pickle import loads, dumps
from six.moves import cStringIO as StringIO
import iotbx.pdb
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long
import os.path
import json
import time

def exercise_cbetadev():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_cbetadev(): input pdb (pdb1jxt.ent) not available")
    return
  from mmtbx.validation import cbetadev
  pdb_in = iotbx.pdb.input(file_name=regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  assert approx_equal(validation.get_weighted_outlier_percent(), 4.40420846587)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    assert (validation.n_outliers == len(validation.results) == 6)
    assert ([ cb.id_str() for cb in validation.results ] ==
      [' A   7 AILE', ' A   8 BVAL', ' A   8 CVAL', ' A  30 BTHR',
       ' A  39 BTHR', ' A  43 BASP'])
    assert approx_equal([ cb.deviation for cb in validation.results ],
      [0.25977096732623106, 0.2577218834868609, 0.6405578498280606,
       0.81238828498566, 0.9239566035292618, 0.5001892640352836])
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    assert not show_diff(out.getvalue(),
"""\
pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:
pdb :A:ile: A:   7 :  0.260: -46.47:   0.45:A:
pdb :B:val: A:   8 :  0.258:  80.92:   0.30:B:
pdb :C:val: A:   8 :  0.641: -53.98:   0.20:C:
pdb :B:thr: A:  30 :  0.812: -76.98:   0.30:B:
pdb :B:thr: A:  39 :  0.924:  56.41:   0.30:B:
pdb :B:asp: A:  43 :  0.500:   7.56:   0.25:B:
SUMMARY: 6 C-beta deviations >= 0.25 Angstrom (Goal: 0)
""")

  # Now with all residues
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    for outlier in validation.results :
      assert (len(outlier.xyz) == 3)
    assert (validation.n_outliers == 6)
    assert (len(validation.results) == 51)
    assert validation.percent_outliers > 0.
    #assert validation.percent_outliers==10.
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    assert not show_diff(out.getvalue(), """\
pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:
pdb : :thr: A:   1 :  0.102:  11.27:   1.00: :
pdb :A:thr: A:   2 :  0.022: -49.31:   0.67:A:
pdb : :cys: A:   3 :  0.038: 103.68:   1.00: :
pdb : :cys: A:   4 :  0.047:-120.73:   1.00: :
pdb : :pro: A:   5 :  0.069:-121.41:   1.00: :
pdb : :ser: A:   6 :  0.052: 112.87:   1.00: :
pdb :A:ile: A:   7 :  0.260: -46.47:   0.45:A:
pdb :B:ile: A:   7 :  0.153: 122.97:   0.55:B:
pdb :A:val: A:   8 :  0.184:-155.36:   0.50:A:
pdb :B:val: A:   8 :  0.258:  80.92:   0.30:B:
pdb :C:val: A:   8 :  0.641: -53.98:   0.20:C:
pdb : :ala: A:   9 :  0.061: -82.84:   1.00: :
pdb :A:arg: A:  10 :  0.023: 172.25:   1.00:A:
pdb : :ser: A:  11 :  0.028:-129.11:   1.00: :
pdb :A:asn: A:  12 :  0.021: -80.80:   0.50:A:
pdb :B:asn: A:  12 :  0.199:  50.01:   0.50:B:
pdb :A:phe: A:  13 :  0.067: -37.32:   0.65:A:
pdb :B:phe: A:  13 :  0.138:  19.24:   0.35:B:
pdb : :asn: A:  14 :  0.065: -96.35:   1.00: :
pdb : :val: A:  15 :  0.138: -96.63:   1.00: :
pdb : :cys: A:  16 :  0.102: -28.64:   1.00: :
pdb : :arg: A:  17 :  0.053:-106.79:   1.00: :
pdb : :leu: A:  18 :  0.053:-141.51:   1.00: :
pdb : :pro: A:  19 :  0.065:-146.95:   1.00: :
pdb : :thr: A:  21 :  0.086:  53.80:   1.00: :
pdb :A:pro: A:  22 :  0.092: -83.39:   0.55:A:
pdb :A:glu: A:  23 :  0.014:-179.53:   0.50:A:
pdb :B:glu: A:  23 :  0.050:-179.78:   0.50:B:
pdb : :ala: A:  24 :  0.056: -88.96:   1.00: :
pdb : :leu: A:  25 :  0.084:-106.42:   1.00: :
pdb : :cys: A:  26 :  0.074: -94.70:   1.00: :
pdb : :ala: A:  27 :  0.056: -62.15:   1.00: :
pdb : :thr: A:  28 :  0.056:-114.82:   1.00: :
pdb :A:tyr: A:  29 :  0.068:   0.22:   0.65:A:
pdb :A:thr: A:  30 :  0.180: 103.27:   0.70:A:
pdb :B:thr: A:  30 :  0.812: -76.98:   0.30:B:
pdb : :cys: A:  32 :  0.029: -84.07:   1.00: :
pdb : :ile: A:  33 :  0.048:-119.17:   1.00: :
pdb : :ile: A:  34 :  0.045:  99.02:   1.00: :
pdb : :ile: A:  35 :  0.052:-128.24:   1.00: :
pdb : :pro: A:  36 :  0.084:-142.29:   1.00: :
pdb : :ala: A:  38 :  0.039:  50.01:   1.00: :
pdb :A:thr: A:  39 :  0.093: -96.63:   0.70:A:
pdb :B:thr: A:  39 :  0.924:  56.41:   0.30:B:
pdb : :cys: A:  40 :  0.013:-144.11:   1.00: :
pdb : :pro: A:  41 :  0.039: -97.09:   1.00: :
pdb :A:asp: A:  43 :  0.130:-146.91:   0.75:A:
pdb :B:asp: A:  43 :  0.500:   7.56:   0.25:B:
pdb : :tyr: A:  44 :  0.085:-143.63:   1.00: :
pdb : :ala: A:  45 :  0.055:  33.32:   1.00: :
pdb : :asn: A:  46 :  0.066: -50.46:   1.00: :
SUMMARY: 6 C-beta deviations >= 0.25 Angstrom (Goal: 0)
""")

  # Auxilary function: extract_atoms_from_residue_group
  from mmtbx.validation.cbetadev import extract_atoms_from_residue_group
  from iotbx import pdb
  pdb_1 = pdb.input(source_info=None, lines="""\
ATOM   1185  N  ASER A 146      24.734  37.097  16.303  0.50 16.64           N
ATOM   1186  N  BSER A 146      24.758  37.100  16.337  0.50 16.79           N
ATOM   1187  CA ASER A 146      24.173  37.500  17.591  0.50 16.63           C
ATOM   1188  CA BSER A 146      24.237  37.427  17.662  0.50 16.87           C
ATOM   1189  C  ASER A 146      22.765  36.938  17.768  0.50 15.77           C
ATOM   1190  C  BSER A 146      22.792  36.945  17.783  0.50 15.94           C
ATOM   1191  O  ASER A 146      22.052  36.688  16.781  0.50 14.91           O
ATOM   1192  O  BSER A 146      22.091  36.741  16.779  0.50 15.17           O
ATOM   1193  CB ASER A 146      24.118  39.035  17.649  0.50 16.93           C
ATOM   1194  CB BSER A 146      24.321  38.940  17.904  0.50 17.48           C
ATOM   1195  OG ASER A 146      23.183  39.485  18.611  0.50 17.56           O
ATOM   1196  OG BSER A 146      23.468  39.645  17.028  0.50 18.32           O  """).construct_hierarchy()
  pdb_2 = pdb.input(source_info=None, lines="""\
ATOM   1185  N   SER A 146      24.734  37.097  16.303  0.50 16.64           N
ATOM   1187  CA  SER A 146      24.173  37.500  17.591  0.50 16.63           C
ATOM   1189  C   SER A 146      22.765  36.938  17.768  0.50 15.77           C
ATOM   1191  O   SER A 146      22.052  36.688  16.781  0.50 14.91           O
ATOM   1193  CB ASER A 146      24.118  39.035  17.649  0.50 16.93           C
ATOM   1194  CB BSER A 146      24.321  38.940  17.904  0.50 17.48           C
ATOM   1195  OG ASER A 146      23.183  39.485  18.611  0.50 17.56           O
ATOM   1196  OG BSER A 146      23.468  39.645  17.028  0.50 18.32           O  """).construct_hierarchy()
  pdb_3 = pdb.input(source_info=None, lines="""\
ATOM   1185  N   SER A 146      24.734  37.097  16.303  0.50 16.64           N
ATOM   1187  CA  SER A 146      24.173  37.500  17.591  0.50 16.63           C
ATOM   1189  C   SER A 146      22.765  36.938  17.768  0.50 15.77           C
ATOM   1191  O   SER A 146      22.052  36.688  16.781  0.50 14.91           O
ATOM   1193  CB  SER A 146      24.118  39.035  17.649  0.50 16.93           C
ATOM   1195  OG ASER A 146      23.183  39.485  18.611  0.50 17.56           O
ATOM   1196  OG BSER A 146      23.468  39.645  17.028  0.50 18.32           O  """).construct_hierarchy()
  rg1 = pdb_1.only_model().only_chain().only_residue_group()
  rg2 = pdb_2.only_model().only_chain().only_residue_group()
  rg3 = pdb_3.only_model().only_chain().only_residue_group()
  all_relevant_atoms_1 = extract_atoms_from_residue_group(rg1)
  all_relevant_atoms_2 = extract_atoms_from_residue_group(rg2)
  all_relevant_atoms_3 = extract_atoms_from_residue_group(rg3)
  keys_1 = [ sorted([ k for k in a.keys() ]) for a in all_relevant_atoms_1 ]
  keys_2 = [ sorted([ k for k in a.keys() ]) for a in all_relevant_atoms_2 ]
  keys_3 = [ sorted([ k for k in a.keys() ]) for a in all_relevant_atoms_3 ]
  assert keys_1 == [[' C  ',' CA ',' CB ',' N  '],[' C  ',' CA ',' CB ',' N  ']]
  assert keys_2 == [[' C  ',' CA ',' CB ',' N  '],[' C  ',' CA ',' CB ',' N  ']]
  assert keys_3 == [[' C  ', ' CA ', ' CB ', ' N  ']]
  print("OK")

def exercise_cbetadev_d_peptide():
  from iotbx import pdb
  from mmtbx.validation import cbetadev
  hierarchy = pdb.input(source_info=None, lines='''
HETATM  566  N   DAL A  11      17.834  32.465  30.842  1.00  5.31           N
HETATM  567  CA  DAL A  11      16.839  31.829  29.896  1.00  4.19           C
HETATM  568  CB  DAL A  11      16.047  30.735  30.590  1.00  5.94           C
HETATM  569  C   DAL A  11      17.623  31.302  28.726  1.00  3.99           C
HETATM  570  O   DAL A  11      18.770  30.854  28.916  1.00  5.73           O''').construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  # assert approx_equal(validation.get_weighted_outlier_percent(), 4.40420846587)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    assert (len(validation.results) == 1)
    assert (validation.n_outliers == 0)
    assert ([ cb.id_str() for cb in validation.results ] ==
      [' A  11  DAL'])
    assert approx_equal([ cb.deviation for cb in validation.results ],
      [0.02848041692354018])
    print(validation.percent_outliers)
    assert validation.percent_outliers == 0.
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    print(out.getvalue())
    for cb in validation.results:
      print(cb.id_str(), cb.deviation)
      assert cb.deviation<1.
  print('OK')

def exercise_cbetadev_misnamed_peptides():
  #testing that residues with wrong chirality show up as outliers
  #I have swapped the names of an otherwise real THR and DTH from 7ooj.pdb
  #they should both be outliers
  from iotbx import pdb
  from mmtbx.validation import cbetadev
  hierarchy = pdb.input(source_info=None, lines='''
HETATM  411  N   THR A  53      22.401  17.450 -18.803  1.00 66.31           N
HETATM  412  CA  THR A  53      23.810  17.544 -18.580  1.00 66.58           C
HETATM  413  CB  THR A  53      24.424  18.930 -18.527  1.00 67.06           C
HETATM  414  CG2 THR A  53      23.940  19.663 -17.276  1.00 68.25           C
HETATM  415  OG1 THR A  53      24.132  19.634 -19.725  1.00 69.19           O
HETATM  416  C   THR A  53      24.317  16.671 -19.722  1.00 67.57           C
HETATM  417  O   THR A  53      25.217  15.882 -19.525  1.00 76.32           O
ATOM    519  N   DTH A  66      14.250  13.187 -35.224  1.00 56.90           N
ATOM    520  CA  DTH A  66      12.879  13.360 -34.677  1.00 58.61           C
ATOM    521  C   DTH A  66      12.864  13.014 -33.188  1.00 57.65           C
ATOM    522  O   DTH A  66      13.234  11.876 -32.846  1.00 58.82           O
ATOM    523  CB  DTH A  66      11.870  12.499 -35.441  1.00 60.31           C
ATOM    524  OG1 DTH A  66      11.977  12.858 -36.818  1.00 66.15           O
ATOM    525  CG2 DTH A  66      10.452  12.683 -34.948  1.00 62.67           C''').construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=True)
  # assert approx_equal(validation.get_weighted_outlier_percent(), 4.40420846587)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    assert (validation.n_outliers == len(validation.results) == 2)
    assert ([ cb.id_str() for cb in validation.results ] ==
      [' A  53  THR',' A  66  DTH'])
    print([ cb.deviation for cb in validation.results ])
    assert approx_equal([ cb.deviation for cb in validation.results ],
      [2.3247665655338596, 2.4166600592687995])
    print(validation.percent_outliers)
    assert validation.percent_outliers>0.
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    print(out.getvalue())
    for cb in validation.results:
      print(cb.id_str(), cb.deviation)
      assert cb.deviation>2.
  print('OK')

def exercise_cbetadev_nonstandard_peptide():
  #testing that a nonstandard animo acid defaults to the general case
  #LYZ is hydroxylysine
  from iotbx import pdb
  from mmtbx.validation import cbetadev
  hierarchy = pdb.input(source_info=None, lines='''
HETATM 3181  N   LYZ D   3      -1.842  -5.028  54.291  1.00 35.06           N
HETATM 3182  CA  LYZ D   3      -3.207  -4.726  53.880  1.00 36.27           C
HETATM 3183  C   LYZ D   3      -3.841  -3.660  54.771  1.00 41.71           C
HETATM 3184  O   LYZ D   3      -4.669  -2.842  54.381  1.00 52.20           O
HETATM 3185  CB  LYZ D   3      -4.176  -5.939  53.836  1.00 34.68           C''').construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  # assert approx_equal(validation.get_weighted_outlier_percent(), 4.40420846587)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    assert (len(validation.results) == 1)
    assert (validation.n_outliers == 0)
    assert ([ cb.id_str() for cb in validation.results ] ==
      [' D   3  LYZ'])
    assert approx_equal([ cb.deviation for cb in validation.results ],
      [0.14107909562037108])
    print(validation.percent_outliers)
    assert validation.percent_outliers == 0.
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    print(out.getvalue())
    for cb in validation.results:
      print(cb.id_str(), cb.deviation)
      assert cb.deviation<1.
  print('OK')

def exercise_cbetadev_unknown_peptide():
  #testing that a nonstandard animo acid defaults to the general case
  #LYZ is hydroxylysine
  from iotbx import pdb
  from mmtbx.validation import cbetadev
  hierarchy = pdb.input(source_info=None, lines='''
HETATM 3181  N   LY? D   3      -1.842  -5.028  54.291  1.00 35.06           N
HETATM 3182  CA  LY? D   3      -3.207  -4.726  53.880  1.00 36.27           C
HETATM 3183  C   LY? D   3      -3.841  -3.660  54.771  1.00 41.71           C
HETATM 3184  O   LY? D   3      -4.669  -2.842  54.381  1.00 52.20           O
HETATM 3185  CB  LY? D   3      -4.176  -5.939  53.836  1.00 34.68           C''').construct_hierarchy()
  validation = cbetadev.cbetadev(
    pdb_hierarchy=hierarchy,
    outliers_only=False)
  # assert approx_equal(validation.get_weighted_outlier_percent(), 4.40420846587)
  for unpickle in [False, True] :
    if unpickle :
      validation = loads(dumps(validation))
    assert (len(validation.results) == 1)
    assert (validation.n_outliers == 0)
    for cb in validation.results: print(cb.id_str())
    assert ([ cb.id_str() for cb in validation.results ] ==
      [' D   3  LY?'])
    assert approx_equal([ cb.deviation for cb in validation.results ],
      [0.14107909562037108])
    print(validation.percent_outliers)
    assert validation.percent_outliers == 0.
    out = StringIO()
    validation.show_old_output(out=out, verbose=True)
    print(out.getvalue())
    for cb in validation.results:
      print(cb.id_str(), cb.deviation)
      assert cb.deviation<1.
  print('OK')

def exercise_cbetadev_json(test_mmcif=False):
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_cbetadev(): input pdb (pdb1jxt.ent) not available")
    return
  from mmtbx.validation import cbetadev
  dm = DataManager()
  if test_mmcif:
    with open(regression_pdb) as f:
      pdb_1jxt_str = f.read()
    pdb_1jxt_str = convert_string_to_cif_long(pdb_1jxt_str, hetatm_name_addition = "", chain_addition="LONGCHAIN")
    dm.process_model_str("1", pdb_1jxt_str)
    m = dm.get_model("1")
    chainA = "ALONGCHAIN"
  else:
    m = dm.get_model(regression_pdb)
    chainA = "A"
  cbeta_json = cbetadev.cbetadev(pdb_hierarchy=m.get_hierarchy(), outliers_only=True).as_JSON()
  cbeta_dict = json.loads(cbeta_json)
  assert len(cbeta_dict['flat_results'])==6, "tst_cbetadev json output not returning correct number of outliers, now: "+str(len(cbeta_dict['flat_results']))
  assert approx_equal(cbeta_dict['flat_results'][0]['deviation'], 0.25977096732623106), "tst_cbetadev json output first deviation not approx_equal, now: "+str(cbeta_dict['flat_results'][0]['deviation'])
  assert approx_equal(cbeta_dict['flat_results'][-1]['deviation'], 0.5001892640352836), "tst_cbetadev json output last deviation not approx_equal, now: "+str(cbeta_dict['flat_results'][-1]['deviation'])
  assert approx_equal(cbeta_dict['hierarchical_results'][''][chainA]["   8 "]['B']["dihedral_NABB"], 80.92016704402938), "tst_cbetadev json output hierarchical result changed dihedral_NABB result, now: "+str(cbeta_dict['hierarchical_results'][''][chainA]["   8 "]['B']["dihedral_NABB"])
  assert cbeta_dict['summary_results'][""]['num_outliers']==6, "tst_cbetadev json output summary results num_outliers changed, now: "+str(cbeta_dict['summary_results'][""]['num_outliers'])
  assert cbeta_dict['summary_results'][""]['num_cbeta_residues']==51, "tst_cbetadev json output summary results num_cbeta_residues changed, now: "+str(cbeta_dict['summary_results'][""]['num_cbeta_residues'])
  return cbeta_dict

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_cbetadev()
  exercise_cbetadev_d_peptide()
  exercise_cbetadev_misnamed_peptides()
  exercise_cbetadev_nonstandard_peptide()
  exercise_cbetadev_unknown_peptide()
  cb_dict = exercise_cbetadev_json()
  cb_dict_cif = exercise_cbetadev_json(test_mmcif=True)
  assert cb_dict['summary_results'] == cb_dict_cif['summary_results'], "tst_cbetadev summary results changed between pdb and cif version"
  print("OK. Time: %8.3f"%(time.time()-t0))
