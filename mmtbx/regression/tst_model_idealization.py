from __future__ import division
from libtbx import easy_run
from mmtbx.secondary_structure.build.tst_2 import tst_01_start_lines
import os.path
import time

def exercise_01(prefix="tst_mi_test_01"):
  # no SS annotations
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(tst_01_start_lines)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      ">& %s.log" % prefix])
  print cmd
  easy_run.call(cmd)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  for l in ["Sorry: No secondary structure annotations found.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()
  # assert os.path.isfile("%s_start.pdb_idealized.pdb" % prefix)

def exercise_02(prefix="tst_mi_test_02"):
  # Same as 01, but with SS annotations in PDB file
  h_records = """\
HELIX    1   1 PRO A    3  ALA A   21  1                                  19
HELIX    2   2 ARG A   23  GLN A   44  1                                  22
"""
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(h_records)
  pdb_file.write(tst_01_start_lines)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      ">& %s.log" % prefix])
  print cmd
  easy_run.call(cmd)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  for l in ["Replacing ss-elements with ideal ones:\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()
  assert os.path.isfile("%s_start.pdb_idealized.pdb" % prefix)
  res_pdb = open("%s_start.pdb_idealized.pdb" % prefix, "r")
  res_pdb_lines = res_pdb.readlines()
  res_pdb.close()
  for l in [
      "HELIX    1   1 PRO A    3  ALA A   21  1                                  19\n",
      "HELIX    2   2 ARG A   23  GLN A   44  1                                  22\n",
      ]:
    assert l in res_pdb_lines, "'%s' not in pdb file." % l
  for l in ["CRYST1", "SCALE1"]:
    found = False
    for pdb_l in res_pdb_lines:
      if pdb_l.startswith(l):
        found = True
        continue
    assert found, "%s not found in pdb file" % l



if (__name__ == "__main__"):
  t0 = time.time()
  exercise_01()
  exercise_02()
  print "Time: %.2f" % (time.time() - t0)
  print "OK"
