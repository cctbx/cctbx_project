from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
import time
from libtbx.utils import null_out

pdb_str = """
CRYST1   20.645   25.843   29.055  90.00  90.00  90.00 P 1
SCALE1      0.048438  0.000000  0.000000        0.00000
SCALE2      0.000000  0.038695  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034417        0.00000
ATOM      1  N   CYS C 284      13.554  13.652   6.082  1.00 43.53      EA   N
ATOM      2  CA  CYS C 284      14.504  14.116   7.088  1.00 43.46      EA   C
ATOM      3  C   CYS C 284      13.857  14.177   8.466  1.00 46.41      EA   C
ATOM      4  O   CYS C 284      14.509  13.899   9.480  1.00 34.73      EA   O
ATOM      5  CB  CYS C 284      15.043  15.493   6.699  1.00 37.37      EA   C
ATOM      6  SG  CYS C 284      15.645  15.628   5.000  1.00 30.55      EA   S
ATOM      7  N   HIS C 285      12.571  14.535   8.520  1.00 81.03      EA   N
ATOM      8  CA  HIS C 285      11.880  14.644   9.799  1.00 78.99      EA   C
ATOM      9  C   HIS C 285      11.736  13.290  10.481  1.00 78.91      EA   C
ATOM     10  O   HIS C 285      11.798  13.206  11.713  1.00 83.81      EA   O
ATOM     11  CB  HIS C 285      10.507  15.286   9.600  1.00 87.56      EA   C
ATOM     12  CG  HIS C 285       9.755  15.507  10.875  1.00100.53      EA   C
ATOM     13  ND1 HIS C 285      10.129  16.454  11.804  1.00102.57      EA   N
ATOM     14  CD2 HIS C 285       8.650  14.905  11.376  1.00110.92      EA   C
ATOM     15  CE1 HIS C 285       9.288  16.426  12.823  1.00110.54      EA   C
ATOM     16  NE2 HIS C 285       8.381  15.495  12.588  1.00115.69      EA   N
ATOM     17  N   GLU C 286      11.547  12.223   9.702  1.00 41.43      EA   N
ATOM     18  CA  GLU C 286      11.389  10.896  10.285  1.00 45.81      EA   C
ATOM     19  C   GLU C 286      12.742  10.245  10.550  1.00 45.31      EA   C
ATOM     20  O   GLU C 286      12.996   9.747  11.652  1.00 47.64      EA   O
ATOM     21  CB  GLU C 286      10.559  10.014   9.350  1.00 31.98      EA   C
ATOM     22  CG  GLU C 286       9.189  10.568   9.001  1.00 31.19      EA   C
ATOM     23  CD  GLU C 286       8.308  10.778  10.211  1.00 46.04      EA   C
ATOM     24  OE1 GLU C 286       7.157  11.227  10.033  1.00 37.93      EA   O
ATOM     25  OE2 GLU C 286       8.766  10.496  11.337  1.00 53.32      EA   O1-
ATOM     26  N   GLY C 289      14.202  12.615  13.228  1.00120.10      EA   N
ATOM     27  CA  GLY C 289      13.414  12.984  14.388  1.00114.20      EA   C
ATOM     28  C   GLY C 289      12.937  11.824  15.240  1.00106.50      EA   C
ATOM     29  O   GLY C 289      12.978  11.900  16.472  1.00135.13      EA   O
ATOM     30  N   HIS C 290      12.487  10.744  14.604  1.00 38.52      EA   N
ATOM     31  CA  HIS C 290      11.863   9.625  15.301  1.00 48.73      EA   C
ATOM     32  C   HIS C 290      12.756   8.398  15.399  1.00 53.95      EA   C
ATOM     33  O   HIS C 290      12.984   7.888  16.500  1.00 51.44      EA   O
ATOM     34  CB  HIS C 290      10.542   9.248  14.614  1.00 62.93      EA   C
ATOM     35  CG  HIS C 290       9.496  10.316  14.686  1.00 56.22      EA   C
ATOM     36  ND1 HIS C 290       9.507  11.425  13.869  1.00 61.73      EA   N
ATOM     37  CD2 HIS C 290       8.406  10.444  15.480  1.00 51.42      EA   C
ATOM     38  CE1 HIS C 290       8.469  12.190  14.155  1.00 65.43      EA   C
ATOM     39  NE2 HIS C 290       7.785  11.618  15.129  1.00 66.69      EA   N
ATOM     40  N   VAL C 291      13.265   7.902  14.267  1.00 73.68      EA   N
ATOM     41  CA  VAL C 291      13.909   6.587  14.249  1.00 64.37      EA   C
ATOM     42  C   VAL C 291      14.982   6.458  15.326  1.00 63.70      EA   C
ATOM     43  O   VAL C 291      14.977   5.453  16.052  1.00 75.32      EA   O
ATOM     44  CB  VAL C 291      14.446   6.282  12.842  1.00 66.63      EA   C
ATOM     45  CG1 VAL C 291      15.263   5.000  12.850  1.00 61.40      EA   C
ATOM     46  CG2 VAL C 291      13.298   6.187  11.848  1.00 83.70      EA   C
ATOM     47  N   VAL C 329       6.501  13.140  22.870  1.00 30.06      EA   N
ATOM     48  CA  VAL C 329       7.335  12.486  21.868  1.00 29.12      EA   C
ATOM     49  C   VAL C 329       7.839  13.498  20.849  1.00 29.42      EA   C
ATOM     50  O   VAL C 329       8.982  13.414  20.387  1.00 30.09      EA   O
ATOM     51  CB  VAL C 329       6.561  11.337  21.192  1.00 26.47      EA   C
ATOM     52  CG1 VAL C 329       7.384  10.727  20.069  1.00 25.53      EA   C
ATOM     53  CG2 VAL C 329       6.187  10.272  22.214  1.00 26.38      EA   C
ATOM     54  N   GLU C 330       7.004  14.476  20.493  1.00 54.98      EA   N
ATOM     55  CA  GLU C 330       7.375  15.437  19.461  1.00 50.88      EA   C
ATOM     56  C   GLU C 330       8.233  16.564  20.026  1.00 53.96      EA   C
ATOM     57  O   GLU C 330       9.287  16.890  19.469  1.00 56.09      EA   O
ATOM     58  CB  GLU C 330       6.115  16.006  18.806  1.00 54.40      EA   C
ATOM     59  CG  GLU C 330       6.312  16.485  17.380  1.00 54.09      EA   C
ATOM     60  CD  GLU C 330       6.616  15.348  16.429  1.00 51.73      EA   C
ATOM     61  OE1 GLU C 330       6.385  14.180  16.807  1.00 45.43      EA   O
ATOM     62  OE2 GLU C 330       7.084  15.620  15.304  1.00 72.06      EA   O1-
ATOM     63  N   PHE C 331       7.801  17.165  21.133  1.00 52.30      EA   N
ATOM     64  CA  PHE C 331       8.480  18.302  21.750  1.00 45.16      EA   C
ATOM     65  C   PHE C 331       8.595  18.100  23.255  1.00 38.54      EA   C
ATOM     66  O   PHE C 331       8.298  18.989  24.055  1.00 40.25      EA   O
ATOM     67  CB  PHE C 331       7.748  19.603  21.437  1.00 51.35      EA   C
ATOM     68  CG  PHE C 331       7.634  19.893  19.970  1.00 50.90      EA   C
ATOM     69  CD1 PHE C 331       8.639  20.575  19.305  1.00 48.59      EA   C
ATOM     70  CD2 PHE C 331       6.523  19.482  19.254  1.00 41.82      EA   C
ATOM     71  CE1 PHE C 331       8.537  20.843  17.953  1.00 56.15      EA   C
ATOM     72  CE2 PHE C 331       6.416  19.747  17.902  1.00 48.23      EA   C
ATOM     73  CZ  PHE C 331       7.423  20.429  17.251  1.00 52.08      EA   C
TER
HETATM   74 FE    FE E   9       6.958  14.164  13.719  1.00129.77          Fe
TER
HETATM   75  O   HOH G  21       5.000  12.310  14.576  1.00 60.96           O
HETATM   76  O   HOH G  22       5.115  15.599  12.863  1.00 38.47           O
HETATM   77  O   HOH G  23       6.308  12.238  11.815  1.00 68.93           O
TER
END"""


def run():
  """
  Check metal linking (make sure CE1 is not included).
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.automatic_linking.link_metals=True
  m = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  m.process(pdb_interpretation_params=params, make_restraints=True)
  grm = m.get_restraints_manager()
  ph = m.get_hierarchy()
  atoms = list(ph.atoms())
  assert atoms[73].element.strip().upper() == "FE"
  bond_proxies_simple, asu = grm.geometry.get_all_bond_proxies(
    sites_cart = ph.atoms().extract_xyz())
  expected_link_atoms = ['NE2', "OE2", "O", "O", "O"]
  expected_link_atoms.sort()
  linked_atoms_found = []
  for p in bond_proxies_simple:
    if(73 in p.i_seqs):
      i,j = p.i_seqs
      a1 = atoms[i].name.strip().upper()
      a2 = atoms[j].name.strip().upper()
      if(i!=73): linked_atoms_found.append(a1)
      else:      linked_atoms_found.append(a2)
  linked_atoms_found.sort()
  assert expected_link_atoms == linked_atoms_found, '%s != %s' % (expected_link_atoms, linked_atoms_found)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.3f"%(time.time()-t0))
