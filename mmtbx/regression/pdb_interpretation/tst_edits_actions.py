from __future__ import absolute_import, division, print_function
import os, time
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

# QUESTION: if we try to change something that is not there shall we Sorry or
# do actualy change "nothing" with something that is provided?
# Depending on the asnwer: check for expected lines using assert_lines_in_file or
# check for expected Sorry in r.stderr_lines[0] .
# ANSWER: For ALL edits the behavior is following:
# add when restraint is there - Sorry
# change when restraint is not there - Sorry.

pdb_str = """\
CRYST1   18.319   14.631   14.035  90.00  90.00  90.00 P 1
ATOM      1  N   ILE A   7       7.541   5.492   8.708  1.00 27.82           N
ATOM      2  CA  ILE A   7       8.118   5.291   7.390  1.00 27.92           C
ATOM      3  C   ILE A   7       9.585   5.691   7.394  1.00 26.95           C
ATOM      4  O   ILE A   7      10.426   5.000   6.823  1.00 27.42           O
ATOM      5  CB  ILE A   7       7.326   6.122   6.344  1.00 28.34           C
ATOM      6  CG1 ILE A   7       5.922   5.540   6.179  1.00 29.61           C
ATOM      7  CG2 ILE A   7       8.046   6.109   5.000  1.00 29.41           C
ATOM      8  CD1 ILE A   7       5.000   6.390   5.343  1.00 30.63           C
TER
END
"""

cmd_base = "phenix.pdb_interpretation write_geo=true %s.pdb %s.eff"

def run_bond(prefix):
  edits = """
refinement {
  geometry_restraints.edits {
    bond {
      action = %s
      atom_selection_1 = %s
      atom_selection_2 = %s
      distance_ideal = 2.0
      sigma = 1.0
    }
  }
}
"""
  sel_old_1 = "chain A and resseq 7 and name C"
  sel_old_2 = "chain A and resseq 7 and name O"
  sel_new_1 = "chain A and resseq 7 and name CA"
  sel_new_2 = "chain A and resseq 7 and name O"
  cntr=0
  for action in ["add","delete","change"]:
    for kind in ["new","old"]:
      prefix_ = "%s_bond_%s_%s"%(prefix, kind, action)
      with open("%s.eff"%prefix_, "w") as f:
        f.write(edits%(action, eval("sel_%s_1"%kind), eval("sel_%s_2"%kind)))
      cmd = cmd_base%(prefix, prefix_)
      print (cmd)
      r = easy_run.fully_buffered(cmd)
      if action == "delete":
        cntr+=1
        assert r.stderr_lines[0] == \
            'Sorry: geometry_restraints.edits.bond.action = delete not implemented.'
      else:
        if kind == "new":
          if action == 'add':
            cntr+=1
            assert_lines_in_file(
                file_name = "%s.pdb.geo"%prefix,
                lines     = """bond pdb=" CA  ILE A   7 "
                             pdb=" O   ILE A   7 "
                          ideal  model  delta    sigma   weight residual
                          2.000  2.394 -0.394 1.00e+00 1.00e+00 1.56e-01""")
          elif action == 'change':
            cntr += 1
            assert r.stderr_lines[0] == "Sorry: Bond below does not exists, use action=add instead."
        elif(kind == "old"):
          if action == 'add':
            cntr += 1
            assert r.stderr_lines[0] == "Sorry: Bond below exists, use action=change instead."
          elif action == 'change':
            cntr+=1
            assert_lines_in_file(
                file_name = "%s.pdb.geo"%prefix,
                lines     = """bond pdb=" C   ILE A   7 "
                                  pdb=" O   ILE A   7 "
                          ideal  model  delta    sigma   weight residual
                          2.000  1.229  0.771 1.00e+00 1.00e+00 5.94e-01""")
  assert cntr==6

def run_angle(prefix):
  edits = """
refinement {
  geometry_restraints.edits {
    angle {
      action = %s
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = 10
      sigma = 1
    }
  }
}
"""
  sel_old_1 = "chain A and resseq 7 and name N"
  sel_old_2 = "chain A and resseq 7 and name CA"
  sel_old_3 = "chain A and resseq 7 and name CB"
  sel_new_1 = "chain A and resseq 7 and name N"
  sel_new_2 = "chain A and resseq 7 and name CA"
  sel_new_3 = "chain A and resseq 7 and name CG2"
  cntr=0
  for action in ["add","delete","change"]:
    for kind in ["new","old"]:
      prefix_ = "%s_angle_%s_%s"%(prefix, kind, action)
      with open("%s.eff"%prefix_, "w") as f:
        f.write(edits%(action, eval("sel_%s_1"%kind),
                               eval("sel_%s_2"%kind),
                               eval("sel_%s_3"%kind)))
      cmd = cmd_base%(prefix, prefix_)
      print (cmd)
      r = easy_run.fully_buffered(cmd)
      if action == "add":
        if kind == "new":
          cntr+=1
          assert_lines_in_file(
              file_name = "%s.pdb.geo"%prefix,
              lines     = """
angle pdb=" N   ILE A   7 "
      pdb=" CA  ILE A   7 "
      pdb=" CG2 ILE A   7 "
    ideal   model   delta    sigma   weight residual
    10.00  143.31 -133.31 1.00e+00 1.00e+00 1.78e+04""")
        elif kind == "old":
          cntr+=1
          assert r.stderr_lines[0] == \
              "Sorry: Some (1) restraints were not added because they are already present."
      if action == "delete":
        cntr+=1
        assert r.stderr_lines[0] == \
            'Sorry: geometry_restraints.edits.angle.action = delete not implemented.'
      if action == "change":
        if kind == "new":
          cntr+=1
          assert r.stderr_lines[0] == \
              'Sorry: Angle below is not restrained, nothing to change.'
        elif kind == "old":
          cntr+=1
          assert_lines_in_file(
              file_name = "%s.pdb.geo"%prefix,
              lines     = """
angle pdb=" N   ILE A   7 "
      pdb=" CA  ILE A   7 "
      pdb=" CB  ILE A   7 "
    ideal   model   delta    sigma   weight residual
    10.00  109.54  -99.54 1.00e+00 1.00e+00 9.91e+03""")
  assert cntr==6


def run(prefix = os.path.basename(__file__).replace(".py","")):
  """
  Test edits: add/delete/change actions for bonds and angles.
  TODO: add dihedral
  """
  with open("%s.pdb"%prefix, "w") as f:
    f.write(pdb_str)
  run_bond(prefix)
  run_angle(prefix)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %5.2f"%(time.time()-t0))
  print("OK")
