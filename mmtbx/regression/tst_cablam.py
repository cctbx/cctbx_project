from __future__ import division
from mmtbx.cablam import cablam_validate
from libtbx.test_utils import show_diff
from iotbx import pdb
import libtbx.load_env
import os

ref_cablam_give_text = """
residue,contour_level,loose_alpha,regular_alpha,loose_beta,regular_beta
ILE    17 ,0.00777,0.00000,0.00000,0.05855,0.00123
ILE    29 ,0.00287,0.00006,0.00000,0.00000,0.00000
ASN    55 ,0.02086,0.00010,0.00000,0.00409,0.00000
PHE   114 ,0.02716,0.00000,0.00000,0.00012,0.00000
LYS   135 ,0.04655,0.00000,0.00000,0.00003,0.00000
"""

ref_cablam_give_oneline = """pdb103l.ent:151:3.3:1.3:0.00
"""

class cablam_test_string():
  #I wrote the regression test to use a class with a custom .write() method as a
  #  proof of principle for learning OOP and to see if I could. Possible because
  #  all my print functions accept an optional writeto= variable.
  def write(self,string):
    self.output += str(string)
  def __init__(self):
    self.output = ""

def exercise_cablam():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb103l.ent",
    test=os.path.isfile) #This is the same file used for tst_kinemage.py
  if (regression_pdb is None):
    print "Skipping exercise_cablam(): input pdb (pdb103l.ent) not available"
    return
  #-----
  pdb_io = pdb.input(regression_pdb)
  pdbid = os.path.basename(regression_pdb)
  hierarchy = pdb_io.construct_hierarchy()

  oneline_test = cablam_test_string()
  cablam_validate.oneline(hierarchy, peptide_cutoff=0.05,
    peptide_bad_cutoff=0.01, ca_cutoff=0.005, pdbid=pdbid,
    writeto=oneline_test)

  text_test = cablam_test_string()
  outliers = cablam_validate.analyze_pdb(hierarchy, outlier_cutoff=0.05, pdbid=pdbid)
  cablam_validate.give_text(outliers, writeto=text_test)

  #print '|',oneline_test.output,'|'
  #print '|',ref_cablam_give_oneline,'|'
  #print '|',text_test.output,'|'
  #print '|',ref_cablam_give_text,'|'

  assert not show_diff(oneline_test.output , ref_cablam_give_oneline)
  assert not show_diff(text_test.output , ref_cablam_give_text)

def run():
  ##if (not libtbx.env.has_module(name="probe")):
  ##  print \
  ##    "Skipping kinemage test:" \
  ##    " probe not available"
  ##elif (not libtbx.env.has_module(name="reduce")):
  ##  print \
  ##    "Skipping kinemage test:" \
  ##    " reduce not available"
  ##else:
  ##  exercise_kinemage()
  ##print "OK"
  exercise_cablam()
  print "OK"

if (__name__ == "__main__"):
  run()
