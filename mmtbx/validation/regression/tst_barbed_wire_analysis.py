from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile

import libtbx.load_env
from libtbx.utils import format_cpu_times
from iotbx.data_manager import DataManager
from mmtbx.validation.barbed_wire_analysis import (
    barbed_wire_analysis, predicted_residue)

CATEGORIES = set(["Predictive", "Near-predictive", "Unpacked high pLDDT",
                  "Pseudostructure", "Barbed wire", "Unphysical"])

# Neither model alone reaches every category. AF_json_v3 is the one with
# "Barbed wire" and "Unphysical", pae_model the one with "Unpacked high pLDDT".
# Both are AlphaFold models with pLDDT in the B-factor column, which is what
# this analysis reads.
MODELS = {
    "AF_json_v3.pdb": {
        "n_residues": 335,
        "expect": set(["Barbed wire", "Unphysical"]),
    },
    "pae_model.pdb": {
        "n_residues": 309,
        "expect": set(["Unpacked high pLDDT"]),
    },
}


def find_model(name):
  return libtbx.env.find_in_repositories(
      relative_path="cctbx_project/mmtbx/regression/pdbs/%s" % name,
      test=os.path.isfile)


class fake_residue(predicted_residue):
  """A residue carrying only what predictalyze() reads, so the decision table
  can be exercised without reduce or probe."""

  def __init__(self, resnum, plddt, packing_quality, signature=False):
    self.chain = 'A'
    self.resseq = '%4d' % resnum
    self.resnum = resnum
    self.resid = 'A,%4d' % resnum
    self.plddt = plddt
    self.packing_quality = packing_quality
    self.barbed_wire_signature = signature
    self.high_outlier_density = None
    self.text_code = ['', '', '', '', '', '']
    self.feedback = ''
    for flag in ('out_rama', 'rama_high_psi', 'rama_high_phi', 'out_omega',
                 'cis_pro', 'out_geom', 'out_geom_cnca', 'out_cablam',
                 'out_ca_geom'):
      setattr(self, flag, None)


def exercise_decision_table():
  """Every branch of predictalyze(), without running the external tools.

  pLDDT 70 and packing quality 1 are the two thresholds; within each half a
  barbed-wire signature or a high outlier density picks the worse of the pair.
  """
  cases = [
      ("Predictive",          85, 2, False),
      ("Unpacked high pLDDT", 85, 0, False),
      ("Near-predictive",     50, 2, False),
      ("Unphysical",          50, 2, True),
      ("Pseudostructure",     50, 0, False),
      ("Barbed wire",         50, 0, True),
  ]
  assert set(c[0] for c in cases) == CATEGORIES

  # Two padding residues at each end, because real chains have them: the
  # 5-residue window in analyze_contacts() sets packing_quality only on indices
  # 2..N-3, and the classify loop skips residues without it, which keeps every
  # classified residue inside categorize_outliers() 3-residue window (1..N-2).
  # Padding here reproduces that, so the cases exercise the decision table.
  residues = [fake_residue(1, 90, 3)]
  for i, (_name, plddt, packing, signature) in enumerate(cases):
    residues.append(fake_residue(i + 2, plddt, packing, signature))
  residues.append(fake_residue(len(cases) + 2, 90, 3))
  residues.append(fake_residue(len(cases) + 3, 90, 3))

  bwa = object.__new__(barbed_wire_analysis)
  bwa.res_list = {'A': residues}
  bwa.res_dict = dict((r.resid, r) for r in residues)
  barbed_wire_analysis.predictalyze(bwa)

  for i, (name, _plddt, _packing, _signature) in enumerate(cases):
    got = residues[i + 1].feedback
    assert got == name, (name, got)

  # 'L' marks low pLDDT and 'p' poor packing, so the codes track the thresholds.
  assert residues[1].text_code == '------', residues[1].text_code
  assert residues[2].text_code == '-p----', residues[2].text_code
  assert residues[3].text_code == 'L-----', residues[3].text_code
  assert residues[6].text_code == 'Lp----', residues[6].text_code


def exercise_model(name, spec):
  """The analysis on a real model.

  Most of what is checked here holds by construction rather than by
  measurement, so it can be asserted exactly without becoming a tripwire for
  another platform or a different probe. The one measured quantity, the share
  of residues carrying contacts, is deliberately a loose floor: it is a
  liveness check on phenix.reduce and phenix.probe, not a claim about how many
  contacts they should find. Exact per-residue expectations belong somewhere
  the environment is pinned, not in this test.
  """
  path = find_model(name)
  if path is None:
    print("Skipping %s: not found" % name)
    return set()

  before = set(os.listdir("."))
  dm = DataManager()
  dm.process_model_file(path)
  bwa = barbed_wire_analysis(model=dm.get_model())
  residues = list(bwa.res_dict.values())
  assert len(residues) == spec["n_residues"], (name, len(residues))

  feedback = [r.feedback for r in residues]
  unknown = set(f for f in feedback if f) - CATEGORIES
  assert not unknown, (name, unknown)
  assert spec["expect"] <= set(feedback), (name, spec["expect"] - set(feedback))

  for chain, chain_residues in bwa.res_list.items():
    n = len(chain_residues)
    classified = [i for i, r in enumerate(chain_residues) if r.feedback]
    scored = [i for i, r in enumerate(chain_residues)
              if r.packing_quality is not None]

    # predictalyze() skips any residue without a packing quality before it
    # assigns anything, so the two sets cannot differ. They did once, when the
    # two sliding windows swept different distances.
    assert classified == scored, (name, chain)

    # analyze_contacts() centres a 5-residue window on res_slice[2], so the
    # first two and last two residues of a chain can never be scored, and
    # nothing else can be missed.
    assert classified == list(range(2, n - 2)), (name, chain, n)

  # Six slots either way, but predictalyze() joins them into a string only for
  # the residues it classifies; the ends keep the list.
  for r in residues:
    assert len(r.text_code) == 6, (name, r.resid, r.text_code)
    expected = str if r.feedback else list
    assert isinstance(r.text_code, expected), (name, r.resid, r.text_code)

  # Loose on purpose. With no contacts at all every packed residue is reported
  # unpacked, which this catches; the exact share is probe's business and moves
  # with its atom typing.
  with_contacts = sum(1 for r in residues
                      if r.packing_hb or r.packing_vdw or r.packing_bo)
  assert with_contacts > len(residues) // 2, (name, with_contacts)

  # do_probe() runs external tools through a file it writes itself.
  left = set(os.listdir(".")) - before
  assert not left, (name, sorted(left))
  return set(f for f in feedback if f)


def main():
  exercise_decision_table()
  cwd = os.getcwd()
  tmp = tempfile.mkdtemp()
  seen = set()
  try:
    os.chdir(tmp)
    for name, spec in sorted(MODELS.items()):
      seen |= exercise_model(name, spec)
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)
  if seen:
    assert seen == CATEGORIES, CATEGORIES - seen
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  main()
