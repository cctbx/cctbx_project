from __future__ import absolute_import, division, print_function
from mmtbx.validation.restraints import chiralities
from libtbx.utils import null_out
from libtbx.easy_pickle import loads
from iotbx.data_manager import DataManager
from mmtbx.model import manager
from libtbx.test_utils import convert_pdb_to_cif_for_pdb_str
import time
import json
import difflib
import io

pdb_2v2k_str = """MODEL        1
ATOM    115  N   DAL A  15      23.144  19.154  24.719  1.00 23.07           N
ATOM    116  CA  DAL A  15      24.113  19.525  23.694  1.00 22.93           C
ATOM    117  C   DAL A  15      23.806  18.844  22.351  1.00 22.65           C
ATOM    118  O   DAL A  15      24.011  19.423  21.279  1.00 22.00           O
ATOM    119  CB  DAL A  15      25.505  19.141  24.158  1.00 22.70           C
HETATM 1607 FE1  F3S A 107      20.751  10.219  14.297  1.00 14.17          FE
HETATM 1608 FE3  F3S A 107      22.460  11.642  15.904  1.00 13.82          FE
HETATM 1609 FE4  F3S A 107      20.669  12.898  14.321  1.00 15.68          FE
HETATM 1610  S1  F3S A 107      22.965   9.822  14.630  1.00 16.10           S
HETATM 1611  S2  F3S A 107      20.493  11.571  12.455  1.00 15.59           S
HETATM 1612  S4  F3S A 107      22.860  13.566  14.695  1.00 16.02           S
HETATM 1613  S3  F3S A 107      20.168  11.468  16.107  1.00 14.13           S
HETATM 1614 FE1  F3S A 108      21.826  14.818  24.986  1.00 15.67          FE
HETATM 1615 FE3  F3S A 108      21.220  13.567  27.239  1.00 15.93          FE
HETATM 1616 FE4  F3S A 108      23.758  13.892  26.554  1.00 16.45          FE
HETATM 1617  S1  F3S A 108      20.202  15.285  26.498  1.00 17.82           S
HETATM 1618  S2  F3S A 108      23.825  15.811  25.405  1.00 16.95           S
HETATM 1619  S3  F3S A 108      22.148  12.719  25.509  1.00 16.89           S
HETATM 1620  S4  F3S A 108      22.801  13.937  28.648  1.00 19.61           S
ENDMDL
MODEL        2
ATOM    115  N   DAL A  15      23.144  19.154  24.719  1.00 23.07           N
ATOM    116  CA  DAL A  15      24.113  19.525  23.694  1.00 22.93           C
ATOM    117  C   DAL A  15      23.806  18.844  22.351  1.00 22.65           C
ATOM    118  O   DAL A  15      24.011  19.423  21.279  1.00 22.00           O
ATOM    119  CB  DAL A  15      25.505  19.141  24.158  1.00 22.70           C
HETATM 1607 FE1  F3S A 107      20.751  10.219  14.297  1.00 14.17          FE
HETATM 1608 FE3  F3S A 107      22.460  11.642  15.904  1.00 13.82          FE
HETATM 1609 FE4  F3S A 107      20.669  12.898  14.321  1.00 15.68          FE
HETATM 1610  S1  F3S A 107      22.965   9.822  14.630  1.00 16.10           S
HETATM 1611  S2  F3S A 107      20.493  11.571  12.455  1.00 15.59           S
HETATM 1612  S4  F3S A 107      22.860  13.566  14.695  1.00 16.02           S
HETATM 1613  S3  F3S A 107      20.168  11.468  16.107  1.00 14.13           S
HETATM 1614 FE1  F3S A 108      21.826  14.818  24.986  1.00 15.67          FE
HETATM 1615 FE3  F3S A 108      21.220  13.567  27.239  1.00 15.93          FE
HETATM 1616 FE4  F3S A 108      23.758  13.892  26.554  1.00 16.45          FE
HETATM 1617  S1  F3S A 108      20.202  15.285  26.498  1.00 17.82           S
HETATM 1618  S2  F3S A 108      23.825  15.811  25.405  1.00 16.95           S
HETATM 1619  S3  F3S A 108      22.148  12.719  25.509  1.00 16.89           S
HETATM 1620  S4  F3S A 108      22.801  13.937  28.648  1.00 19.61           S
ENDMDL
END
"""

expected_chiral = """  atoms                   ideal    model    delta   sigma  residual   deviation
   A 107  F3S  S2
   A 107  F3S FE1
   A 107  F3S FE3
   A 107  F3S FE4         10.77     7.78     2.99  2.00e-01  2.23e+02  14.9*sigma
   A 107  F3S  S2
   A 107  F3S FE1
   A 107  F3S FE3
   A 107  F3S FE4         10.77     7.78     2.99  2.00e-01  2.23e+02  14.9*sigma
   A  15  DAL  CA
   A  15  DAL  N
   A  15  DAL  C
   A  15  DAL  CB         -2.48     2.54    -5.02  2.00e-01  6.31e+02  25.1*sigma
   A  15  DAL  CA
   A  15  DAL  N
   A  15  DAL  C
   A  15  DAL  CB         -2.48     2.54    -5.02  2.00e-01  6.31e+02  25.1*sigma
   A 108  F3S  S2
   A 108  F3S FE1
   A 108  F3S FE3
   A 108  F3S FE4         10.77    -6.76    17.53  2.00e-01  7.68e+03  87.7*sigma
   A 108  F3S  S2
   A 108  F3S FE1
   A 108  F3S FE3
   A 108  F3S FE4         10.77    -6.76    17.53  2.00e-01  7.68e+03  87.7*sigma

  Min. delta:    2.988
  Max. delta:   17.532
  Mean delta:   10.670
"""

def calculate_results():
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("1",pdb_2v2k_str)
  model = dm.get_model("1")
  model.set_stop_for_unknowns(False)
  p = manager.get_default_pdb_interpretation_params()
  ##print(dir(p.pdb_interpretation))
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.flip_symmetric_amino_acids=False
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
  model.set_log(log = null_out())
  model.process(make_restraints=True, pdb_interpretation_params=p)
  geometry_restraints_manager = model.get_restraints_manager().geometry
  pdb_hierarchy = model.get_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xray_structure = model.get_xray_structure()
  from mmtbx import restraints
  restraints_manager = restraints.manager(
    geometry=geometry_restraints_manager)
  sites_cart = xray_structure.sites_cart()
  hd_selection = xray_structure.hd_selection()
  pdb_atoms = pdb_hierarchy.atoms()
  energies_sites = restraints_manager.energies_sites(
    sites_cart=sites_cart,
    compute_gradients=False).geometry
  restraint_proxies = getattr(restraints_manager.geometry, "chirality_proxies")
  chiral_list = chiralities(
      pdb_atoms=pdb_atoms,
      sites_cart=sites_cart,
      energies_sites=energies_sites,
      restraint_proxies=restraint_proxies,
      unit_cell=xray_structure.unit_cell(),
      ignore_hd=True,
      sigma_cutoff=4.0,
      outliers_only=True,
      use_segids_in_place_of_chainids=False)
  return chiral_list

def exercise_chiral_validation(chiral_list):

  chiral_io = io.StringIO()
  chiral_list.show(out=chiral_io, verbose=True)
  chiral_results = chiral_io.getvalue()
  chiral_io.close()
  diff = difflib.unified_diff(chiral_results.splitlines(), expected_chiral.splitlines(), fromfile="testvalue", tofile="expectedvalue")
  changed_lines = ""
  for line in diff:
    if line.startswith("-") or line.startswith("+"):
      changed_lines = changed_lines+"\n"+line
  assert changed_lines == "", "chiral validation text output changed, at the following lines: "+changed_lines

def exercise_chiral_json(chiral_list):
  ch_dict = json.loads(chiral_list.as_JSON())
  #import pprint
  #pprint.pprint(ch_dict)
  assert len(ch_dict['flat_results']) == 6, "tst_chiral_validation json output not returning correct number of water clashes, now: "+str(len(ch_dict['flat_results']))
  assert ch_dict['flat_results'][0]["outlier_type"] == "Tetrahedral geometry outlier", "tst_chiral_validation json output first outlier_type value changed, now: "+ch_dict['flat_results'][0]["outlier_type"]
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(ch_dict['hierarchical_results'], "Tetrahedral geometry outlier")==2, "tst_chiral_validation json hierarchical output total number of Tetrahedral geometry outlier changed, now: "+str(test_utils.count_dict_values(ch_dict['hierarchical_results'], "Tetrahedral geometry outlier"))
  assert test_utils.count_dict_values(ch_dict['hierarchical_results'], "Pseudochiral naming error")==2, "tst_chiral_validation json hierarchical output total number of Pseudochiral naming error changed, now: "+str(test_utils.count_dict_values(ch_dict['hierarchical_results'], "Pseudochiral naming error"))
  assert test_utils.count_dict_values(ch_dict['hierarchical_results'], "Chiral handedness swap")==2, "tst_chiral_validation json hierarchical output total number of Chiral handedness swap changed, now: "+str(test_utils.count_dict_values(ch_dict['hierarchical_results'], "Chiral handedness swap"))
  summary_results_dict = ch_dict['summary_results']
  if "   1" in summary_results_dict:
    summary_results_1_dict = summary_results_dict["   1"]
    summary_results_2_dict = summary_results_dict["   2"]
  else:
    summary_results_1_dict = summary_results_dict["1"]
    summary_results_2_dict = summary_results_dict["2"]
  assert summary_results_1_dict["num_outliers"] == 3, "tst_chiral_validation json summary output total number of outliers changed, now: "+str(ch_dict['summary_results']["   1"]["num_outliers"])
  assert summary_results_1_dict["num_chiral_centers"] == 1, "tst_chiral_validation json summary output total number of true chiral centers changed, now: "+str(ch_dict['summary_results']["   1"]["num_chiral_centers"])
  assert summary_results_1_dict["num_total"] == 3, "tst_chiral_validation json summary output total number of tetrahedral changed, now: "+str(ch_dict['summary_results']["   1"]["num_total"])
  assert summary_results_2_dict["num_outliers"] == 3, "tst_chiral_validation json summary output model 2 total number of outliers changed, now: "+str(ch_dict['summary_results']["   1"]["num_outliers"])
  return ch_dict

if (__name__ == "__main__"):
  t0 = time.time()
  chiral_list = calculate_results()
  exercise_chiral_validation(chiral_list)
  ch_dict = exercise_chiral_json(chiral_list)
  convert_pdb_to_cif_for_pdb_str(locals(), chain_addition="LONGCHAIN", hetatm_name_addition = "", key_str="pdb_", print_new_string = False)
  chiral_list = calculate_results()
  ch_dict_cif = exercise_chiral_json(chiral_list)
  assert ch_dict['summary_results']['   1'] == ch_dict_cif['summary_results']['1'], "tst_chiral_validation summary results changed between pdb and cif version"
  assert ch_dict['summary_results']['   2'] == ch_dict_cif['summary_results']['2'], "tst_chiral_validation summary results changed between pdb and cif version"
  print("OK. Time: %8.3f"%(time.time()-t0))
