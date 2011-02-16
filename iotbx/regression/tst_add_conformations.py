
from libtbx.test_utils import contains_lines, Exception_expected
from libtbx.utils import Sorry
import libtbx.load_env
import cStringIO
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not available, skipping test."
    return
  from iotbx.command_line.pdb_add_conformations import run
  out = cStringIO.StringIO()
  run([pdb_file], out=out)
  assert contains_lines(out.getvalue(), "Modified model: 4254 atoms")
  out = cStringIO.StringIO()
  run([pdb_file, "atom_selection=\"chain A and not resname HOH\""], out=out)
  assert contains_lines(out.getvalue(), "Modified model: 3990 atoms")
  run([pdb_file, "new_occ=0.4", "atom_selection=\"resseq 1:275\""], out=out)
  from iotbx import file_reader
  pdb_in = file_reader.any_file("1ywf_split.pdb", force_type="pdb").file_object
  atoms = pdb_in.atoms()
  occ = atoms.extract_occ()
  assert (occ.count(0.6) == occ.count(0.4) == 1858)
  out = cStringIO.StringIO()
  run([pdb_file, "n_confs=3", "new_occ=0.25"], out=out)
  pdb_in = file_reader.any_file("1ywf_split.pdb", force_type="pdb").file_object
  assert contains_lines(out.getvalue(), """\
WARNING: zero-occupancy atom:
HETATM 1941  O  AHOH A 354      -0.009  56.525  -3.872  0.25 29.17           O\
""")
  atoms = pdb_in.atoms()
  assert (atoms.size() == 6381)
  occ = atoms.extract_occ()
  assert (occ.count(0.5) == 2126) and (occ.count(0.25) == 4254)
  try :
    run([pdb_file, "atom_selection=\"chain G\""], out=out)
  except Sorry, e :
    assert (str(e) == "Empty selection.")
  else :
    raise Exception_expected
  try :
    run([pdb_file, "new_occ=2"], out=out)
  except Sorry, e :
    assert (str(e) == "new_occ must be between 0 and 1.0")
  else :
    raise Exception_expected
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1akg.pdb",
    test=os.path.isfile)
  try :
    run([pdb_file], out=out)
  except Sorry, e :
    assert (str(e) == """\
Atom group included in selection already has one or more alternate conformers:
ATOM     22  OG ASER     4      -1.752   0.849   3.272  0.50 11.67           O\
""")
  else :
    raise Exception_expected
  run([pdb_file, "atom_selection=\"not name OG\""], out=out)
  print "OK"

if (__name__ == "__main__") :
  exercise()
