from __future__ import absolute_import, division, print_function
from mmtbx.validation import ramalyze
from iotbx.data_manager import DataManager
import time

# ------------------------------------------------------------------------------
# Testing for edge cases in Ramalyze
#  Put more edge cases here if convenient
# ------------------------------------------------------------------------------

#Twisted prolines (omega >30 degrees from planar) should be evaluated as cisPro or TransPro
#  depending on whether omega is closer to cis or trans
#Abundant twisted Pro are likely to be encountered in low-pLDDT AF2 predictions
pdb_str_1 = """
ATOM     88  N   LEU A  13     -24.165 -36.507 -21.753  1.00 55.67           N
ATOM     89  CA  LEU A  13     -24.423 -36.560 -23.221  1.00 55.67           C
ATOM     90  C   LEU A  13     -23.595 -35.472 -23.965  1.00 55.67           C
ATOM     91  CB  LEU A  13     -25.921 -36.323 -23.558  1.00 55.67           C
ATOM     92  O   LEU A  13     -23.431 -34.383 -23.423  1.00 55.67           O
ATOM     93  CG  LEU A  13     -26.799 -37.590 -23.681  1.00 55.67           C
ATOM     94  CD1 LEU A  13     -28.274 -37.188 -23.721  1.00 55.67           C
ATOM     95  CD2 LEU A  13     -26.491 -38.387 -24.949  1.00 55.67           C
ATOM     96  N   PRO A  14     -23.223 -35.648 -25.252  1.00 50.20           N
ATOM     97  CA  PRO A  14     -22.490 -36.757 -25.871  1.00 50.20           C
ATOM     98  C   PRO A  14     -21.129 -36.303 -26.478  1.00 50.20           C
ATOM     99  CB  PRO A  14     -23.437 -37.243 -26.971  1.00 50.20           C
ATOM    100  O   PRO A  14     -20.833 -35.114 -26.574  1.00 50.20           O
ATOM    101  CG  PRO A  14     -24.150 -35.964 -27.431  1.00 50.20           C
ATOM    102  CD  PRO A  14     -23.914 -34.937 -26.316  1.00 50.20           C
ATOM    103  N   GLY A  15     -20.286 -37.252 -26.906  1.00 53.60           N
ATOM    104  CA  GLY A  15     -18.938 -36.976 -27.444  1.00 53.60           C
ATOM    105  C   GLY A  15     -18.873 -36.479 -28.909  1.00 53.60           C
ATOM    106  O   GLY A  15     -19.866 -36.559 -29.634  1.00 53.60           O
ATOM    516  N   SER A  72      24.352  38.834 -48.364  1.00 40.52           N
ATOM    517  CA  SER A  72      25.319  39.477 -47.445  1.00 40.52           C
ATOM    518  C   SER A  72      26.747  39.550 -48.012  1.00 40.52           C
ATOM    519  CB  SER A  72      25.429  38.689 -46.125  1.00 40.52           C
ATOM    520  O   SER A  72      27.107  38.729 -48.855  1.00 40.52           O
ATOM    521  OG  SER A  72      24.180  38.209 -45.657  1.00 40.52           O
ATOM    522  N   PRO A  73      27.603  40.420 -47.444  1.00 48.63           N
ATOM    523  CA  PRO A  73      29.011  40.085 -47.181  1.00 48.63           C
ATOM    524  C   PRO A  73      29.415  40.292 -45.700  1.00 48.63           C
ATOM    525  CB  PRO A  73      29.835  40.954 -48.142  1.00 48.63           C
ATOM    526  O   PRO A  73      28.603  40.711 -44.878  1.00 48.63           O
ATOM    527  CG  PRO A  73      28.848  41.982 -48.704  1.00 48.63           C
ATOM    528  CD  PRO A  73      27.586  41.810 -47.858  1.00 48.63           C
ATOM    529  N   GLY A  74      30.658  39.927 -45.353  1.00 40.95           N
ATOM    530  CA  GLY A  74      31.201  39.925 -43.981  1.00 40.95           C
ATOM    531  C   GLY A  74      32.034  41.164 -43.582  1.00 40.95           C
ATOM    532  O   GLY A  74      32.067  42.142 -44.325  1.00 40.95           O
"""

def exercise_ramalyze():
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("1",pdb_str_1)
  model = dm.get_model("1")
  hierarchy = model.get_hierarchy()
  r = ramalyze.ramalyze(hierarchy)
  pro14 = r.results[0]
  pro73 = r.results[1]
  assert(pro14.residue_type() == "cis-Pro") #omega= -54.24
  assert(pro73.residue_type() == "trans-Pro") #omega= 138.82

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_ramalyze()
  print("OK. Time: %8.3f"%(time.time()-t0))
