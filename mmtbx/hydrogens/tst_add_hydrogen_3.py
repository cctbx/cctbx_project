from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out
#from libtbx.test_utils import approx_equal

# ------------------------------------------------------------------------------

def run():
  test_000(pdb_str = pdb_str_000)
  test_001(pdb_str = pdb_str_001)

# ------------------------------------------------------------------------------

def test_000(pdb_str):
  '''
    Make sure reduce does not crash for single_atom_residue models
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  # initial model (has no H atoms)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  number_h_expected = model_initial.get_hd_selection().count(True)
  assert(number_h_expected == 0)
  # place H atoms
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_initial)
  reduce_add_h_obj.run()
  # We don't expect H atoms to be placed
  # (not enough restraints for single atom residues)
  model_h_added = reduce_add_h_obj.get_model()
  number_h_placed = model_h_added.get_hd_selection().count(True)
  assert(number_h_placed == 0)

# ------------------------------------------------------------------------------

def test_001(pdb_str):
  '''
    Check keyword n_terminal_charge:
    NH3 on resseq 1, first residue in chain, or no NH3 at all
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  # initial model
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  #
  # place H atoms: NH3 at resseq 1 only
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_initial)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==1)
  assert(h_names_added.count(' H2 ')==1)
  assert(h_names_added.count(' H3 ')==1)
  #
  # place H atoms: NH3 at first residue in chain
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    n_terminal_charge = 'first_in_chain')
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==3)
  assert(h_names_added.count(' H2 ')==3)
  assert(h_names_added.count(' H3 ')==3)
  #
  # place H atoms: no NH3
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    n_terminal_charge = 'no_charge')
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==0)
  assert(h_names_added.count(' H2 ')==0)
  assert(h_names_added.count(' H3 ')==0)

# ------------------------------------------------------------------------------

pdb_str_000 = """
REMARK Make sure reduce does not crash for single_atom_residue models
CRYST1   22.029   33.502   24.035  90.00  90.00  90.00 P 1
ATOM      6  P     A A  10     -62.272  56.445  13.820  1.00 15.00           P
ATOM      7  P     G A  11     -63.673  51.410  11.026  1.00 15.00           P
ATOM      8  P     U A  12     -62.888  45.926   9.711  1.00 15.00           P
ATOM      9  P     U A  13     -60.326  41.305  11.244  1.00 15.00           P
ATOM     10  P     U A  14     -57.909  36.481  13.207  1.00 15.00           P
ATOM     11  P     G A  15     -62.106  32.943  15.800  1.00 15.00           P
ATOM     12  P     A A  16     -65.446  37.240  15.291  1.00 15.00           P
ATOM     13  P     U A  17     -66.286  42.354  18.232  1.00 15.00           P
ATOM     14  P     C A  18     -64.629  46.517  21.258  1.00 15.00           P
ATOM     15  P     A A  19     -60.460  50.019  23.746  1.00 15.00           P
ATOM     16  P     U A  20     -54.257  51.133  23.481  1.00 15.00           P
"""

pdb_str_001 = """
REMARK Make sure NH3 is applied correctly depending on keyword
CRYST1   30.200   47.800   61.300  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   TYR A   5       8.831  48.837  54.788  1.00  9.67           N
ATOM      2  CA  TYR A   5       7.706  49.436  54.084  1.00  9.24           C
ATOM      3  C   TYR A   5       6.456  48.981  54.822  1.00 10.02           C
ATOM      4  O   TYR A   5       6.310  47.784  55.139  1.00 10.62           O
ATOM      5  CB  TYR A   5       7.599  48.942  52.633  1.00  9.83           C
ATOM      6  CG  TYR A   5       8.692  49.472  51.736  1.00 12.79           C
ATOM      7  CD1 TYR A   5       9.959  48.898  51.781  1.00 13.60           C
ATOM      8  CD2 TYR A   5       8.415  50.537  50.880  1.00 12.13           C
ATOM      9  CE1 TYR A   5      10.961  49.400  50.960  1.00 14.77           C
ATOM     10  CE2 TYR A   5       9.426  51.032  50.065  1.00 14.21           C
ATOM     11  CZ  TYR A   5      10.685  50.467  50.116  1.00 14.05           C
ATOM     12  OH  TYR A   5      11.708  50.978  49.318  1.00 17.48           O
ATOM     13  H1  TYR A   5       9.284  48.311  54.231  1.00  9.67           H
ATOM     14  H2  TYR A   5       9.367  49.479  55.092  1.00  9.67           H
ATOM     15  H3  TYR A   5       8.530  48.354  55.472  1.00  9.67           H
ATOM     16  HA  TYR A   5       7.805  50.400  54.049  1.00  9.24           H
ATOM     17  HB2 TYR A   5       7.653  47.974  52.627  1.00  9.83           H
ATOM     18  HB3 TYR A   5       6.748  49.229  52.266  1.00  9.83           H
ATOM     19  HD1 TYR A   5      10.133  48.187  52.354  1.00 13.60           H
ATOM     20  HD2 TYR A   5       7.564  50.912  50.855  1.00 12.13           H
ATOM     21  HE1 TYR A   5      11.811  49.024  50.976  1.00 14.77           H
ATOM     22  HE2 TYR A   5       9.255  51.741  49.488  1.00 14.21           H
ATOM     23  HH  TYR A   5      12.418  50.550  49.452  1.00 17.48           H
TER
ATOM     24  N   GLU B   1      14.684  52.510  58.829  1.00 14.47           N
ATOM     25  CA  GLU B   1      13.950  53.593  58.225  1.00 15.28           C
ATOM     26  C   GLU B   1      12.566  53.009  57.937  1.00 14.74           C
ATOM     27  O   GLU B   1      12.459  51.916  57.371  1.00 14.27           O
ATOM     28  CB  GLU B   1      14.667  53.949  56.967  1.00 18.89           C
ATOM     29  CG  GLU B   1      13.973  54.945  56.083  1.00 27.57           C
ATOM     30  CD  GLU B   1      14.729  55.326  54.802  1.00 32.66           C
ATOM     31  OE1 GLU B   1      15.802  54.776  54.508  1.00 34.50           O
ATOM     32  OE2 GLU B   1      14.224  56.201  54.090  1.00 36.48           O
ATOM     33  H1  GLU B   1      15.398  52.322  58.332  1.00 14.47           H
ATOM     34  H2  GLU B   1      14.945  52.748  59.646  1.00 14.47           H
ATOM     35  H3  GLU B   1      14.162  51.791  58.882  1.00 14.47           H
ATOM     36  HA  GLU B   1      13.870  54.399  58.759  1.00 15.28           H
ATOM     37  HB2 GLU B   1      15.529  54.326  57.204  1.00 18.89           H
ATOM     38  HB3 GLU B   1      14.790  53.139  56.447  1.00 18.89           H
ATOM     39  HG2 GLU B   1      13.835  55.760  56.590  1.00 27.57           H
ATOM     40  HG3 GLU B   1      13.119  54.573  55.814  1.00 27.57           H
TER
ATOM     41  N   PHE C  -3       6.984  40.342  58.778  1.00  8.80           N
ATOM     42  CA  PHE C  -3       7.384  40.247  60.166  1.00  8.32           C
ATOM     43  C   PHE C  -3       7.719  38.788  60.513  1.00  9.15           C
ATOM     44  O   PHE C  -3       8.710  38.555  61.201  1.00  9.31           O
ATOM     45  CB  PHE C  -3       6.284  40.754  61.091  1.00  8.94           C
ATOM     46  CG  PHE C  -3       6.705  40.642  62.560  1.00  9.84           C
ATOM     47  CD1 PHE C  -3       7.825  41.288  63.022  1.00 10.62           C
ATOM     48  CD2 PHE C  -3       5.989  39.828  63.426  1.00 12.63           C
ATOM     49  CE1 PHE C  -3       8.229  41.132  64.328  1.00 11.20           C
ATOM     50  CE2 PHE C  -3       6.398  39.679  64.737  1.00 13.74           C
ATOM     51  CZ  PHE C  -3       7.527  40.326  65.197  1.00 12.55           C
ATOM     52  H1  PHE C  -3       6.155  40.662  58.729  1.00  8.80           H
ATOM     53  H2  PHE C  -3       7.538  40.888  58.346  1.00  8.80           H
ATOM     54  H3  PHE C  -3       7.013  39.534  58.405  1.00  8.80           H
ATOM     55  HA  PHE C  -3       8.167  40.801  60.311  1.00  8.32           H
ATOM     56  HB2 PHE C  -3       6.101  41.686  60.894  1.00  8.94           H
ATOM     57  HB3 PHE C  -3       5.483  40.224  60.960  1.00  8.94           H
ATOM     58  HD1 PHE C  -3       8.313  41.834  62.449  1.00 10.62           H
ATOM     59  HD2 PHE C  -3       5.232  39.381  63.123  1.00 12.63           H
ATOM     60  HE1 PHE C  -3       8.988  41.578  64.629  1.00 11.20           H
ATOM     61  HE2 PHE C  -3       5.909  39.138  65.314  1.00 13.74           H
ATOM     62  HZ  PHE C  -3       7.808  40.220  66.077  1.00 12.55           H
TER
END
"""

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
