from __future__ import absolute_import, division, print_function
import time
from libtbx import easy_pickle
import libtbx.load_env
import mmtbx.idealized_aa_residues.rotamer_manager as rm

def exercise():
  """
  Rotamer chi angle tables are loaded on demand, once per residue type.
  """
  path = libtbx.env.find_in_repositories("chem_data/rotamer_chi_angles")
  o = rm.load(rotamers = "favored")
  assert len(o.result) == 0 # nothing is read until asked for
  ser = o.get_chi_angles(resname = " SER ")
  assert list(o.result.keys()) == ["ser"]
  assert ser == easy_pickle.load("%s/ser_favored.pkl"%path)
  assert o.get_chi_angles(resname = "ser") is ser # read once
  # preload: names that are not amino acids are ignored
  o.preload(resnames = ["THR", "HOH", "ser", "thr"])
  assert sorted(o.result.keys()) == ["ser", "thr"]
  # explicit list of residues: read at construction, as before
  o = rm.load(rotamers = "favored_allowed", residues = ["VAL"])
  assert list(o.result.keys()) == ["val"]
  assert o.get_chi_angles("val") == \
    easy_pickle.load("%s/val_favored_allowed.pkl"%path)
  # no table for this one (file does not exist): None, as before
  assert o.get_chi_angles("ala") is None
  # survives pickling (rotamer manager is passed to parallel workers)
  o.rotamer_evaluator.log = None # sys.stdout cannot be pickled
  o2 = easy_pickle.loads(easy_pickle.dumps(o))
  assert sorted(o2.result.keys()) == ["ala", "val"]
  assert o2.get_chi_angles("gly") is None

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
