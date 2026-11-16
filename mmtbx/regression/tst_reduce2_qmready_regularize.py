"""reduce2 with qmready=True regularizes the geometry of the finished model
(mmtbx.model.manager.regularize_geometry: bond outliers and clashes, with
the smallest shift) once the H are placed and optimized, before the clash
check and before the caps. Models: 1yjp Asn 3 and Gln 4 with C3-N4 = 3.1 A,
so qmready adds OXT to Asn 3 and NH3 to Gln 4 and OXT ASN A 3 - H3 GLN A 4 is
0.892 A (a clash the check used to stop on); and a single Ala, which gets
its NH3 and OXT without a violation. Both in a P 1 box without symmetry
contacts.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.refinement.geometry_minimization import geometry_violations

pdb_str_gap = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM     13  N   ASN A   3      23.342  28.731  29.139  1.00 12.26           N
ATOM     14  CA  ASN A   3      24.587  29.045  29.823  1.00 11.74           C
ATOM     15  C   ASN A   3      25.825  28.473  29.129  1.00 11.10           C
ATOM     16  O   ASN A   3      25.908  27.260  28.882  1.00 10.42           O
ATOM     17  CB  ASN A   3      24.521  28.519  31.276  1.00 12.15           C
ATOM     18  CG  ASN A   3      25.774  28.880  32.095  1.00 12.82           C
ATOM     19  OD1 ASN A   3      26.078  30.066  32.306  1.00 15.05           O
ATOM     20  ND2 ASN A   3      26.509  27.856  32.540  1.00 13.48           N
ATOM     21  N   GLN A   4      28.025  30.547  28.442  1.00 10.29           N
ATOM     22  CA  GLN A   4      29.414  30.207  28.043  1.00 10.53           C
ATOM     23  C   GLN A   4      30.465  30.925  28.932  1.00 10.24           C
ATOM     24  O   GLN A   4      30.577  32.162  28.959  1.00  8.86           O
ATOM     25  CB  GLN A   4      29.686  30.467  26.555  1.00  9.80           C
ATOM     26  CG  GLN A   4      30.974  29.777  26.057  1.00 10.25           C
ATOM     27  CD  GLN A   4      31.534  30.363  24.755  1.00 12.43           C
ATOM     28  OE1 GLN A   4      31.774  31.587  24.654  1.00 14.62           O
ATOM     29  NE2 GLN A   4      31.780  29.480  23.753  1.00  9.05           N
END
"""

pdb_str_clean = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      20.000  20.000  20.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      21.458  20.000  20.000  1.00 10.00           C
ATOM      3  C   ALA A   1      21.985  21.418  20.000  1.00 10.00           C
ATOM      4  O   ALA A   1      21.197  22.361  20.000  1.00 10.00           O
ATOM      5  CB  ALA A   1      21.968  19.248  18.785  1.00 10.00           C
END
"""

def run(model_text, prefix, extra_args=()):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  fn_desc = prefix + "_desc.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in [fn_out, fn_desc]:
    if os.path.exists(fn): os.remove(fn)
  args = [fn_in, "output.filename=" + fn_out,
    "output.description_file_name=" + fn_desc, "overwrite=True"]
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=args + list(extra_args))
  assert os.path.isfile(fn_out), fn_out
  assert os.path.isfile(fn_desc), fn_desc
  return fn_out, open(fn_desc).read()

def get_model(file_name):
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(file_name=file_name), log=null_out())
  # the X-H distances of qmready
  p = model.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = True
  model.process(make_restraints=True, pdb_interpretation_params=p)
  return model

def violations(file_name):
  model = get_model(file_name)
  v = geometry_violations(
    geometry_restraints_manager = model.get_restraints_manager().geometry,
    sites_cart                  = model.get_sites_cart())
  return v.bonds + v.clashes

def heavy_atoms(text):
  """id_str -> xyz of the non-H atoms."""
  result = {}
  for a in iotbx.pdb.input(source_info=None, lines=text).atoms():
    if a.element.strip().upper() in ["H", "D"]: continue
    result[a.id_str()] = a.xyz
  return result

def max_heavy_shift(text_1, text_2):
  """Largest distance between the common non-H atoms of two models."""
  h1, h2 = heavy_atoms(text_1), heavy_atoms(text_2)
  common = set(h1) & set(h2)
  assert len(common) > 0
  result = 0
  for k in common:
    d = sum((a-b)**2 for a, b in zip(h1[k], h2[k]))**0.5
    result = max(result, d)
  return result

def test_gap_regularized():
  """The clash the check used to stop on is regularized away: the model is
  written, free of bond outliers and clashes, with a small shift."""
  fn_out, desc = run(pdb_str_gap, "tst_reduce2_qmready_regularize_1",
    ["qmready=True"])
  assert len(violations(fn_out)) == 0, violations(fn_out)
  assert "Regularize geometry" in desc, desc
  shift = max_heavy_shift(pdb_str_gap, open(fn_out).read())
  assert 0.01 < shift <= 1.0, shift
  # The regularized model needs nothing on a second run.
  fn_out_2, desc_2 = run(open(fn_out).read(),
    "tst_reduce2_qmready_regularize_1b", ["qmready=True"])
  assert "regularization not needed" in desc_2, desc_2
  assert max_heavy_shift(open(fn_out).read(), open(fn_out_2).read()) == 0

def test_clean_untouched():
  """No violations: the non-H atoms are the input ones, and the description
  says so."""
  fn_out, desc = run(pdb_str_clean, "tst_reduce2_qmready_regularize_2",
    ["qmready=True"])
  assert len(violations(fn_out)) == 0, violations(fn_out)
  assert "regularization not needed" in desc, desc
  assert max_heavy_shift(pdb_str_clean, open(fn_out).read()) == 0

def test_no_qmready():
  """Without qmready nothing is regularized."""
  fn_out, desc = run(pdb_str_gap, "tst_reduce2_qmready_regularize_3")
  assert "Regularize geometry" not in desc, desc
  assert "regularization not needed" not in desc, desc
  assert max_heavy_shift(pdb_str_gap, open(fn_out).read()) == 0

def test_capping():
  """With capping=True the regularized coordinates are the capped ones."""
  fn_out, desc = run(pdb_str_gap, "tst_reduce2_qmready_regularize_4",
    ["qmready=True", "capping=True"])
  assert "Regularize geometry" in desc, desc
  fn_ref = "tst_reduce2_qmready_regularize_1_out.pdb"
  assert max_heavy_shift(open(fn_ref).read(), open(fn_out).read()) == 0

if __name__ == "__main__":
  test_gap_regularized()
  test_clean_untouched()
  test_no_qmready()
  test_capping()
  print("OK")
