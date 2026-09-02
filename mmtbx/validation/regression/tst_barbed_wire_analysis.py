from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile

import libtbx.load_env
from libtbx.utils import format_cpu_times
from iotbx.data_manager import DataManager
from mmtbx.validation import barbed_wire_analysis

# pae_model.pdb is a real AlphaFold model, 309 residues with pLDDT in the
# B-factor column, which is what this analysis is for.
MODEL = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/regression/pdbs/pae_model.pdb",
    test=os.path.isfile)

CATEGORIES = set(["Predictive", "Near-predictive", "Unpacked high pLDDT",
                  "Pseudostructure", "Barbed wire"])

# do_probe() shells out to phenix.reduce and phenix.probe through a file it
# writes into the current directory under a fixed name, so the test runs
# somewhere private and checks the run leaves nothing behind.


def run_analysis():
  dm = DataManager()
  dm.process_model_file(MODEL)
  return barbed_wire_analysis.barbed_wire_analysis(model=dm.get_model())


def exercise():
  before = set(os.listdir("."))
  bwa = run_analysis()
  residues = list(bwa.res_dict.values())
  assert len(residues) == 309, len(residues)

  feedback = [r.feedback for r in residues]
  unknown = set(f for f in feedback if f) - CATEGORIES
  assert not unknown, unknown

  # The two residues at each end have no classification: the backbone measures
  # this depends on need neighbours on both sides.
  blank = sorted(int(r.resseq) for r in residues if not r.feedback)
  assert blank == [1, 2, 307, 308, 309], blank

  # A confident model should come out mostly predictive. This is the assertion
  # that fails loudly if phenix.probe or phenix.reduce cannot be run: with no
  # contacts every packed residue is reported unpacked instead.
  predictive = feedback.count("Predictive")
  assert predictive >= 250, predictive

  # Same failure, checked directly rather than through the classification.
  with_contacts = sum(1 for r in residues
                      if r.packing_hb or r.packing_vdw or r.packing_bo)
  assert with_contacts >= 200, with_contacts

  for name in ("Near-predictive", "Unpacked high pLDDT"):
    assert name in feedback, name

  # Six slots, one per problem class, '-' where that problem is absent.
  for r in residues:
    code = r.text_code if isinstance(r.text_code, str) else "".join(r.text_code)
    assert len(code) in (0, 6), (r.resid, repr(code))

  left = set(os.listdir(".")) - before
  assert not left, "files left behind in %s: %s" % (os.getcwd(), sorted(left))


def main():
  if MODEL is None:
    print("Skipping: pae_model.pdb not available")
    return
  cwd = os.getcwd()
  tmp = tempfile.mkdtemp()
  try:
    os.chdir(tmp)
    exercise()
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  main()
