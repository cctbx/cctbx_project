
from __future__ import division
from __future__ import print_function
from libtbx.test_utils import approx_equal, show_diff
import libtbx.load_env
from libtbx.easy_pickle import loads, dumps
from cStringIO import StringIO
import os.path

def exercise_cbetadev():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_cbetadev(): input pdb (pdb1jxt.ent) not available")
    return
  from mmtbx.validation import cbetadev
  from iotbx import file_reader
  pdb_in = file_reader.any_file(file_name=regression_pdb)
  hierarchy = pdb_in.file_object.hierarchy
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

if (__name__ == "__main__"):
  exercise_cbetadev()
