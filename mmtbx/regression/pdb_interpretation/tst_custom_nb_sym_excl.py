from __future__ import absolute_import, division, print_function

from libtbx.test_utils import assert_lines_in_file, assert_lines_not_in_file, assert_lines_in_text
from libtbx import easy_run
import libtbx.load_env
import os

def exercise_01(prefix="tst_custom_nb_sym_excl"):
  fname = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/regression/pdbs/1yjp_h.pdb",
    test=os.path.isfile)

  import mmtbx
  from pathlib import Path
  data_dir = Path(mmtbx.__file__).parent / 'regression' / 'pdbs'
  fname = str( data_dir / '1yjp_h.pdb')

  with open(fname, 'r') as f:
    f_content = f.read()
    f_content = f_content.replace("3.905  1.00 12.26", "3.905  0.33 12.26")
    with open("%s.pdb" % prefix, 'w') as f2:
      f2.write(f_content)

  critical='  nonbonded pdb=" OD1 ASN A   3 "'

  cmd = "mmtbx.pdb_interpretation %s.pdb custom_nonbonded_symmetry_exclusions='resname HOH' custom_nonbonded_symmetry_exclusions='resname HOH and resseq 9' " % prefix
  fb = easy_run.fully_buffered(cmd)
  assert not fb.return_code
  stdout_text = "\n".join(fb.stdout_lines)
  assert_lines_in_text(stdout_text, """\
      WARNING: 1 atom in previous symmetry exclusion group
        Example: pdb=" O   HOH A   9 "
      WARNING: 1 atom with full occupancy
        Example: pdb=" O   HOH A   9 " """)
  assert_lines_in_text(stdout_text, critical)

  cmd = "mmtbx.pdb_interpretation %s.pdb custom_nonbonded_symmetry_exclusions='chain A and resid 3' > %s_2.log" % (prefix, prefix)
  fb = easy_run.fully_buffered(cmd)
  assert not fb.return_code

  assert_lines_not_in_file("%s_2.log" % prefix, critical)
  assert_lines_in_file("%s_2.log" % prefix, """\
      Custom nonbonded symmetry exclusions:
        Atom selection: chain A and resid 3
          Number of atoms selected: 14
          WARNING: 13 atoms with full occupancy
            Example: pdb=" CA  ASN A   3 "
          """)

if(__name__ == "__main__"):
  if libtbx.env.find_in_repositories(relative_path="chem_data") is None:
    print("Skipping exercise_01(): chem_data directory not available")
  else:
    exercise_01()
    print('OK')
