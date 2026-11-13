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

# The second zinc site of 1nlt as deposited (residues 159-162 and 185-188 of
# chain A, ZN 352): Zn-S 1.42-1.45 A (ideal 2.33 A), so the S of the four Cys
# are 2.3-2.5 A from each other. In a P 1 box.
pdb_str_zn = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM    378  N   CYS A 159      12.013  13.368  10.570  1.00110.29           N
ATOM    379  CA  CYS A 159      12.488  13.897  11.844  1.00113.35           C
ATOM    380  C   CYS A 159      11.702  15.147  12.219  1.00113.94           C
ATOM    381  O   CYS A 159      11.781  16.169  11.540  1.00114.95           O
ATOM    382  CB  CYS A 159      13.986  14.218  11.752  1.00116.26           C
ATOM    383  SG  CYS A 159      14.814  14.596  13.336  1.00122.23           S
ATOM    384  N   THR A 160      10.941  15.051  13.304  1.00114.88           N
ATOM    385  CA  THR A 160      10.127  16.161  13.786  1.00114.84           C
ATOM    386  C   THR A 160      10.985  17.322  14.281  1.00114.30           C
ATOM    387  O   THR A 160      10.692  18.489  14.013  1.00113.64           O
ATOM    388  CB  THR A 160       9.207  15.702  14.934  1.00115.91           C
ATOM    389  OG1 THR A 160       9.984  15.007  15.919  1.00116.25           O
ATOM    390  CG2 THR A 160       8.118  14.777  14.409  1.00116.35           C
ATOM    391  N   SER A 161      12.048  16.982  15.003  1.00113.50           N
ATOM    392  CA  SER A 161      12.964  17.967  15.557  1.00112.98           C
ATOM    393  C   SER A 161      13.846  18.593  14.483  1.00112.07           C
ATOM    394  O   SER A 161      14.596  19.531  14.752  1.00111.32           O
ATOM    395  CB  SER A 161      13.837  17.305  16.621  1.00112.41           C
ATOM    396  OG  SER A 161      13.038  16.685  17.612  1.00113.39           O
ATOM    397  N   CYS A 162      13.752  18.065  13.267  1.00111.36           N
ATOM    398  CA  CYS A 162      14.535  18.562  12.141  1.00110.27           C
ATOM    399  C   CYS A 162      13.791  18.246  10.844  1.00106.35           C
ATOM    400  O   CYS A 162      14.004  17.199  10.234  1.00106.26           O
ATOM    401  CB  CYS A 162      15.922  17.902  12.137  1.00114.48           C
ATOM    402  SG  CYS A 162      15.946  16.075  11.960  1.00121.48           S
ATOM    581  N   CYS A 185      20.109  15.969  15.296  1.00126.58           N
ATOM    582  CA  CYS A 185      18.769  16.437  15.605  1.00128.99           C
ATOM    583  C   CYS A 185      18.093  15.342  16.425  1.00131.38           C
ATOM    584  O   CYS A 185      17.803  14.266  15.912  1.00133.18           O
ATOM    585  CB  CYS A 185      17.986  16.695  14.312  1.00128.16           C
ATOM    586  SG  CYS A 185      16.235  16.142  14.294  1.00127.20           S
ATOM    587  N   ASP A 186      17.869  15.615  17.706  1.00134.03           N
ATOM    588  CA  ASP A 186      17.228  14.661  18.612  1.00134.62           C
ATOM    589  C   ASP A 186      17.419  13.200  18.219  1.00135.17           C
ATOM    590  O   ASP A 186      18.411  12.578  18.576  1.00135.49           O
ATOM    591  CB  ASP A 186      15.729  14.966  18.717  1.00134.62           C
ATOM    592  CG  ASP A 186      15.398  15.874  19.882  1.00135.80           C
ATOM    593  OD1 ASP A 186      15.990  16.967  19.979  1.00135.45           O
ATOM    594  OD2 ASP A 186      14.536  15.492  20.700  1.00135.27           O
ATOM    595  N   VAL A 187      16.453  12.666  17.482  1.00135.07           N
ATOM    596  CA  VAL A 187      16.475  11.281  17.032  1.00134.83           C
ATOM    597  C   VAL A 187      16.945  11.143  15.583  1.00134.02           C
ATOM    598  O   VAL A 187      16.464  10.275  14.848  1.00132.66           O
ATOM    599  CB  VAL A 187      15.063  10.642  17.165  1.00135.31           C
ATOM    600  CG1 VAL A 187      15.170   9.124  17.103  1.00136.47           C
ATOM    601  CG2 VAL A 187      14.411  11.097  18.458  1.00136.35           C
ATOM    602  N   CYS A 188      17.870  12.001  15.168  1.00131.92           N
ATOM    603  CA  CYS A 188      18.374  11.947  13.807  1.00130.02           C
ATOM    604  C   CYS A 188      19.812  12.481  13.724  1.00129.02           C
ATOM    605  O   CYS A 188      20.504  12.583  14.740  1.00128.75           O
ATOM    606  CB  CYS A 188      17.436  12.725  12.855  1.00129.20           C
ATOM    607  SG  CYS A 188      17.277  14.534  13.113  1.00126.91           S
TER
HETATM 1815 ZN    ZN A 352      16.052  15.301  13.150  1.00124.81          ZN
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
  with open(fn_desc) as fh:
    desc = fh.read()
  return fn_out, desc

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
  with open(fn_out) as fh:
    out_text = fh.read()
  shift = max_heavy_shift(pdb_str_gap, out_text)
  assert 0.01 < shift <= 1.0, shift
  # The regularized model needs nothing on a second run.
  fn_out_2, desc_2 = run(out_text,
    "tst_reduce2_qmready_regularize_1b", ["qmready=True"])
  assert "regularization not needed" in desc_2, desc_2
  with open(fn_out_2) as fh:
    out_text_2 = fh.read()
  assert max_heavy_shift(out_text, out_text_2) == 0

def test_clean_untouched():
  """No violations: the non-H atoms are the input ones, and the description
  says so."""
  fn_out, desc = run(pdb_str_clean, "tst_reduce2_qmready_regularize_2",
    ["qmready=True"])
  assert len(violations(fn_out)) == 0, violations(fn_out)
  assert "regularization not needed" in desc, desc
  with open(fn_out) as fh:
    assert max_heavy_shift(pdb_str_clean, fh.read()) == 0

def test_no_qmready():
  """Without qmready nothing is regularized."""
  fn_out, desc = run(pdb_str_gap, "tst_reduce2_qmready_regularize_3")
  assert "Regularize geometry" not in desc, desc
  assert "regularization not needed" not in desc, desc
  with open(fn_out) as fh:
    assert max_heavy_shift(pdb_str_gap, fh.read()) == 0

def test_capping():
  """With capping=True the regularized coordinates are the capped ones."""
  fn_out, desc = run(pdb_str_gap, "tst_reduce2_qmready_regularize_4",
    ["qmready=True", "capping=True"])
  assert "Regularize geometry" in desc, desc
  fn_ref = "tst_reduce2_qmready_regularize_1_out.pdb"
  with open(fn_ref) as fh:
    ref_text = fh.read()
  with open(fn_out) as fh:
    out_text = fh.read()
  assert max_heavy_shift(ref_text, out_text) == 0

def test_metal_site():
  """Collapsed bonds to a metal (1nlt: Zn-S 1.42-1.45 A) are bond outliers
  too: qmready regularizes them to within 0.2 A of the ideal 2.33 A, the S
  stay thiolates and nothing needs a cap."""
  fn_out, desc = run(pdb_str_zn, "tst_reduce2_qmready_regularize_5",
    ["qmready=True"])
  assert "Regularize geometry" in desc, desc
  assert len(violations(fn_out)) == 0, violations(fn_out)
  hierarchy = iotbx.pdb.input(file_name=fn_out).construct_hierarchy()
  atoms = hierarchy.atoms()
  zn = [a for a in atoms if a.element.strip().upper() == "ZN"]
  assert len(zn) == 1
  sg = [a for a in atoms if a.name.strip() == "SG"]
  assert len(sg) == 4
  for a in sg:
    assert 2.13 < a.distance(zn[0]) < 2.5, (a.id_str(), a.distance(zn[0]))
  for i, a in enumerate(sg):
    for b in sg[i+1:]:
      assert a.distance(b) > 3.2, (a.id_str(), b.id_str(), a.distance(b))
  for a in sg:
    names = [b.name.strip() for b in a.parent().atoms()]
    assert "HG" not in names, a.id_str()
  assert "Capped 0 atoms" in desc, desc

if __name__ == "__main__":
  test_gap_regularized()
  test_clean_untouched()
  test_no_qmready()
  test_capping()
  test_metal_site()
  print("OK")
