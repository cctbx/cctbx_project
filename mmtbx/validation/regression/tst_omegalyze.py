from __future__ import absolute_import, division, print_function
from mmtbx.validation import omegalyze
from libtbx.test_utils import show_diff
from iotbx import pdb
import libtbx.load_env
import os

ref_omegalyze_give_text = """residues:type:omega:conformation:mc_bmax
 A  40  PHE to  A  41  PRO: Pro     : -14.27:Cis     :28.76
 A 206  VAL to  A 207  LEU: non-Pro : 123.05:Twisted :29.35
 A 504  ILE to  A 505  PRO: Pro     :  -0.31:Cis     :27.55
 A 584  ASN to  A 585  LYS: non-Pro : -12.68:Cis     :26.39
 A 603  ILE to  A 604  GLY: non-Pro :-144.84:Twisted :27.51
 B 929  ARG to  B 930  LEU: non-Pro :-136.14:Twisted :26.97
 B1331  LYS to  B1332  ASP: non-Pro : 134.92:Twisted :27.81
 B1474  GLU to  B1475  LYS: non-Pro : -19.25:Cis     :43.95
 B1593  LYS to  B1594  PRO: Pro     :   8.04:Cis     :26.47
SUMMARY: 3 cis prolines out of 77 PRO
SUMMARY: 0 twisted prolines out of 77 PRO
SUMMARY: 2 other cis residues out of 1464 nonPRO
SUMMARY: 4 other twisted residues out of 1464 nonPRO
"""

class omegalyze_test_string():
  #I wrote the regression test to use a class with a custom .write() method as a
  #  proof of principle for learning OOP and to see if I could. Possible because
  #  all my print functions accept an optional writeto= variable.
  def write(self,string):
    self.output += str(string)
  def __init__(self):
    self.output = ""

def exercise_omegalyze():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2hr0.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_omegalyze(): input pdb (2hr0.pdb) not available")
    return
  #-----
  pdb_io = pdb.input(regression_pdb)
  pdbid = os.path.basename(regression_pdb)
  hierarchy = pdb_io.construct_hierarchy()

  text_test = omegalyze_test_string()
  outliers = omegalyze.omegalyze(
    pdb_hierarchy=hierarchy,
    nontrans_only=True,
    out=text_test,
    quiet=False)
  outliers.show_old_output(out=text_test, verbose=True)

  assert not show_diff(text_test.output , ref_omegalyze_give_text)

def run():
  exercise_omegalyze()
  print("OK")

if (__name__ == "__main__"):
  run()
