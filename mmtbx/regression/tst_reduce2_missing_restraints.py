"""
reduce2 and residues without restraints: the Sorry of the default run (the
residues, why, and what to do), and the log and description wording with
ignore_missing_restraints=True. The sites are cut from 1cc1 (FCO), 1c4c (HC0)
and 4u9i (EJ2): metal cofactors whose CCD entry RDKit does not accept.
"""
from __future__ import absolute_import, division, print_function
import os
import re
import subprocess
import sys
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
import mmtbx.chemical_components
from mmtbx.chemical_components import get_cif_filename
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens import tst_add_hydrogen_19 as t19
from mmtbx.regression import reduce2_neutron_utils as u

FIRST = u.MISSING_RESTRAINTS

# 1cc1: Cys L 73 and the Fe(CN)2(CO) of the [NiFeSe] site; no H in the CCD
pdb_fco = """\
CRYST1  110.390   63.700   99.580  90.00  90.00  90.00 P 21 21 21
ATOM   2536  N   CYS L  73      40.878  11.931  43.670  1.00 19.73           N
ATOM   2537  CA  CYS L  73      40.630  10.641  43.027  1.00 18.38           C
ATOM   2538  C   CYS L  73      39.535  10.695  41.950  1.00 19.80           C
ATOM   2539  O   CYS L  73      38.676   9.799  41.879  1.00 19.76           O
ATOM   2540  CB  CYS L  73      40.196   9.648  44.087  1.00 18.07           C
ATOM   2541  SG  CYS L  73      38.549  10.008  44.757  1.00 22.54           S
HETATM 5890 FE   FCO L 499      37.915   8.881  46.642  1.00 21.64          FE
HETATM 5891  C1  FCO L 499      37.605   8.324  48.404  1.00 22.28           C
HETATM 5892  N1  FCO L 499      37.561   7.902  49.504  1.00 16.70           N
HETATM 5893  C2  FCO L 499      39.589   8.016  46.622  1.00 23.66           C
HETATM 5894  N2  FCO L 499      40.563   7.357  46.576  1.00 20.67           N
HETATM 5895  C3  FCO L 499      36.799   7.586  45.825  1.00 28.98           C
HETATM 5896  O3  FCO L 499      36.116   6.653  45.546  1.00 28.98           O
END
"""

# 1c4c: Cys A 503 and the 2Fe subcluster; the CCD has 2 H (a water ligand)
pdb_hc0 = """\
CRYST1  133.500   83.800   55.600  90.00  90.00  90.00 P 21 21 2
ATOM   3876  N   CYS A 503      53.994  25.962  14.532  1.00 12.04           N
ATOM   3877  CA  CYS A 503      54.500  24.623  14.741  1.00 12.91           C
ATOM   3878  C   CYS A 503      55.099  23.976  13.502  1.00 13.55           C
ATOM   3879  O   CYS A 503      55.987  23.129  13.615  1.00 11.52           O
ATOM   3880  CB  CYS A 503      53.381  23.763  15.273  1.00 12.47           C
ATOM   3881  SG  CYS A 503      52.889  24.099  16.992  1.00 13.37           S
HETATM 4463 FE1  HC0 A 580      51.683  22.203  17.717  1.00 17.44          FE
HETATM 4464 FE2  HC0 A 580      49.549  21.931  19.173  1.00 19.36          FE
HETATM 4465  S1  HC0 A 580      49.707  23.348  17.350  1.00 16.18           S
HETATM 4466  S2  HC0 A 580      51.311  23.213  19.876  1.00 19.37           S
HETATM 4467  O1  HC0 A 580      50.076  25.405  19.007  1.00 17.05           O
HETATM 4468  O8  HC0 A 580      47.780  23.943  20.339  1.00 29.09           O
HETATM 4469  O3  HC0 A 580      54.172  20.853  18.528  1.00 17.75           O
HETATM 4470  O4  HC0 A 580      51.592  21.050  14.981  1.00 17.39           O
HETATM 4471  O5  HC0 A 580      50.612  19.395  17.923  1.00 19.94           O
HETATM 4472  O6  HC0 A 580      49.636  20.238  21.621  1.00 17.50           O
HETATM 4473  O7  HC0 A 580      47.296  20.401  18.084  1.00 11.14           O
HETATM 4474  C3  HC0 A 580      53.233  21.372  18.209  1.00 15.20           C
HETATM 4475  C4  HC0 A 580      51.649  21.479  16.035  1.00 15.17           C
HETATM 4476  C5  HC0 A 580      50.795  20.399  18.390  1.00 17.08           C
HETATM 4477  C6  HC0 A 580      49.578  20.903  20.706  1.00 18.44           C
HETATM 4478  C7  HC0 A 580      48.140  20.994  18.493  1.00 14.98           C
HETATM 4479  C8  HC0 A 580      48.465  23.169  19.909  1.00 21.78           C
END
"""

# 4u9i: Cys L 81 and the [NiFe] site with its deposited hydride (H)
pdb_ej2 = """\
CRYST1   63.850   99.690  122.190  90.00  90.00  90.00 P 21 21 21
ATOM   5166  N   CYS L  81      13.947   9.108 -14.326  1.00  5.29           N
ATOM   5167  CA  CYS L  81      14.109   7.962 -15.205  1.00  5.23           C
ATOM   5168  C   CYS L  81      12.955   6.982 -15.027  1.00  5.55           C
ATOM   5169  O   CYS L  81      12.646   6.577 -13.902  1.00  6.20           O
ATOM   5170  CB  CYS L  81      15.441   7.266 -14.889  1.00  5.75           C
ATOM   5171  SG  CYS L  81      15.843   5.951 -16.054  1.00  5.47           S
HETATM12808 FE   EJ2 L 601      19.249   5.075 -13.055  1.00  5.26          FE2+
HETATM12809 NI   EJ2 L 601      17.788   5.315 -15.164  1.00  5.66          NI2+
HETATM12810  C1  EJ2 L 601      21.026   5.409 -13.514  1.00  5.75           C
HETATM12811  N1  EJ2 L 601      22.114   5.592 -13.838  1.00  6.31           N
HETATM12812  C2  EJ2 L 601      19.620   3.270 -12.628  1.00  5.66           C
HETATM12813  N2  EJ2 L 601      19.817   2.145 -12.312  1.00  5.84           N
HETATM12814  C3  EJ2 L 601      19.499   5.427 -11.383  1.00  5.33           C
HETATM12815  O3  EJ2 L 601      19.663   5.636 -10.236  1.00  6.68           O
HETATM12816  H   EJ2 L 601      19.079   4.553 -14.681  1.00  3.00           H0
END
"""

cryst1 = "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"

def has_ccd(*names):
  for name in names:
    if not os.path.isfile(get_cif_filename(name)): return False
  return True

def sorry_text(prefix, pdb_str, **params):
  try:
    u.run_prog(prefix, pdb_str, **params)
  except Sorry as e:
    return str(e)
  raise AssertionError("no Sorry: %s" % prefix)

def copies_04(n):
  # n residues 04 (no dictionary anywhere), 5 A apart
  lines = []
  for i in range(n):
    for l in t19.pdb_04.splitlines():
      lines.append("%s%4d    %8.3f%s" % (l[:22], i + 1,
        float(l[30:38]) + 5 * i, l[38:]))
  return cryst1.replace("30.000   30.000   30.000",
    "90.000   30.000   30.000") + "\n".join(lines) + "\nEND\n"

def run_flag(prefix, pdb_str, **params):
  """-> (log, description, stdout, H/D names by resname) of a run with the
  flag"""
  out, sys.stdout = sys.stdout, StringIO()
  try:
    prog = u.run_prog(prefix, pdb_str, ignore_missing_restraints=True,
      **params)
    printed = sys.stdout.getvalue()
  finally:
    sys.stdout = out
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  hd = {}
  with open(prog.params.output.filename) as f:
    for l in f:
      if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() in ("H", "D"):
        hd.setdefault(l[17:20].strip(), []).append(l[12:16].strip())
  return prog.log_text, description, printed, hd

class edited_ccd(object):
  """with edited_ccd(name, edit): the CCD entry of name is edit(its text)"""
  def __init__(self, name, edit):
    self.name, self.edit = name, edit
  def __enter__(self):
    cc = mmtbx.chemical_components
    with open(get_cif_filename(self.name)) as f:
      text = self.edit(f.read())
    path = os.path.abspath("tst_r2_missing_%s.cif" % self.name)
    with open(path, "w") as f:
      f.write(text)
    self.original = cc.get_cif_filename
    cc.get_cif_filename = lambda code: path if code.strip() == self.name \
      else self.original(code)
    cc.loaded_cifs.pop(self.name, None)
  def __exit__(self, *args):
    cc = mmtbx.chemical_components
    cc.get_cif_filename = self.original
    cc.loaded_cifs.pop(self.name, None)

def no_type_column(text):
  assert "_chem_comp_atom.type_symbol" in text
  return text.replace("_chem_comp_atom.type_symbol", "_chem_comp_atom.type_x")

def no_atom_loop(text):
  result, n = re.subn(r"loop_\n_chem_comp_atom\.comp_id.*?\n#\n", "", text,
    flags=re.S)
  assert n == 1
  return result

def unknown_h_type(mark):
  def edit(text):
    result, n = re.subn(r"(HC0 +HO1\d +HO1\d +)H ", r"\g<1>%s " % mark, text)
    assert n == 2
    return result
  return edit

def blank_element(pdb_str, name):
  # the element and charge columns of the atom lines of name: blank
  lines = []
  for l in pdb_str.splitlines():
    if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == name:
      l = l[:76]
    lines.append(l)
  return "\n".join(lines) + "\n"

def two_models(pdb_str, ids=(7, 42)):
  lines = pdb_str.splitlines()
  atoms = "".join(l + "\n" for l in lines if l.startswith(("ATOM", "HETATM")))
  return lines[0] + "\n" + "".join("MODEL     %4d\n%sENDMDL\n" % (i, atoms)
    for i in ids) + "END\n"

def to_unl(pdb_str):
  return pdb_str.replace(" 04 B", "UNL B")
# ------------------------------------------------------------------------------

def test_stop_ccd_entry_not_usable():
  if not has_ccd("FCO"):
    print("Skipping test_stop_ccd_entry_not_usable: no CCD"); return
  text = sorry_text("tst_r2_missing_1", pdb_fco)
  lines = text.splitlines()
  assert lines[0] == FIRST + " FCO", text
  assert "  FCO L 499" in lines, text
  assert ("  FCO: the CCD entry exists but could not be used (Explicit "
    "valence for atom # 6 O, 3, is greater than permitted)") in lines, text
  assert text.count("Explicit valence") == 1, text
  assert "Supply a restraint CIF" in text, text
  for part in ("ignore_missing_restraints=True", "gets no H",
               "keep_existing_H=True", "no bonds and no hydrogen-bond types",
               "scored as clashes"):
    assert part in text, (part, text)

def test_stop_no_entry():
  text = sorry_text("tst_r2_missing_2", cryst1 + t19.pdb_04)
  lines = text.splitlines()
  assert lines[0] == FIRST + " 04", text
  assert "  04 B 1" in lines, text
  assert ("  04: no entry in GeoStd, the monomer library or the CCD, and no "
    "restraint file for it") in lines, text
  assert "could not be used" not in text, text

def test_stop_optimize_says_unchanged():
  # approach=optimize: the flag leaves the residue as in the input
  text = sorry_text("tst_r2_missing_3", cryst1 + t19.pdb_04,
    approach="optimize")
  assert text.splitlines()[0] == FIRST + " 04", text
  assert "left exactly as in the input" in text, text
  assert "keep_existing_H" not in text, text

def test_stop_long_list_is_capped():
  text = sorry_text("tst_r2_missing_4", copies_04(13))
  lines = text.splitlines()
  assert lines[0] == FIRST + " 04", text
  listed = [l for l in lines if l.startswith("  04 B ")]
  assert listed == ["  04 B %d" % i for i in range(1, 11)], text
  assert "  ... and 3 more of 04 (13 in all)" in lines, text

def test_stop_names_model_and_altloc():
  # two models; the residue shares its residue group with a water in B
  ag = t19.pdb_04.replace("  04 B", "A 04 B") + (
    "HETATM   60  O  BHOH B   1      15.000  12.000  15.000  0.50 20.00"
    "           O\n")
  # the MODEL IDs are the input's
  text = sorry_text("tst_r2_missing_5", two_models(cryst1 + ag))
  lines = text.splitlines()
  assert "  04 B 1 altloc A (model 7)" in lines, text
  assert "  04 B 1 altloc A (model 42)" in lines, text

def test_stop_unl():
  # a restraint CIF is not used for UNL (mon_lib_query())
  text = sorry_text("tst_r2_missing_11", cryst1 + to_unl(t19.pdb_04))
  lines = text.splitlines()
  assert lines[0] == FIRST + " UNL", text
  assert "Supply a restraint CIF" not in text, text
  assert ("  UNL: a restraint CIF for UNL is not used. Identify the component "
    "and give it its") in lines, text
  assert ("  component ID in the model and in the restraint file."
    ) in lines, text
  assert "ignore_missing_restraints=True" in text, text
  # with another name: the CIF advice is for that one
  both = copies_04(2).splitlines()
  both = "\n".join([to_unl(l) if l[22:26] == "   2" else l for l in both])
  text = sorry_text("tst_r2_missing_12", both + "\n")
  lines = text.splitlines()
  assert lines[0] == FIRST + " 04 UNL", text
  assert "  Supply a restraint CIF for 04." in lines, text
  assert "  UNL: a restraint CIF for UNL is not used. Identify the " \
    "component and give it its" in lines, text

def test_flag_unl_advice():
  # the description gives UNL the advice of the stop; another residue
  # without restraints keeps the restraint CIF
  cys = "".join(l + "\n" for l in pdb_fco.splitlines() if " CYS " in l)
  unl = "UNL: a restraint CIF for UNL is not used. Identify the component " \
    "and give it its component ID in the model and in the restraint file."
  cif = "  What to do: Supply a restraint CIF for the residue."
  log, description, printed, hd = run_flag("tst_r2_missing_19",
    cryst1 + cys + to_unl(t19.pdb_04))
  assert "  What to do: " + unl in description.splitlines(), description
  assert "Supply a restraint CIF" not in description, description
  both = copies_04(2).splitlines()
  both = [to_unl(l) if l[22:26] == "   2" else l for l in both]
  log, description, printed, hd = run_flag("tst_r2_missing_20",
    "\n".join(both[:-1]) + "\n" + cys + "END\n")
  lines = description.splitlines()
  i = lines.index(cif)
  assert lines[i + 1] == "  " + unl, description
  assert "  04 B 1: all H (" in description, description
  assert "  UNL B 2: all H (" in description, description

def test_no_rdkit_needed():
  # a model whose restraints all come from the library, in a fresh process
  # where rdkit cannot be imported
  code = """
import sys
sys.path[:0] = %r
sys.modules["rdkit"] = None
from mmtbx.regression import reduce2_neutron_utils as u
prog = u.run_prog("tst_r2_missing_21", u.pdb_1yjp)
with open(prog.params.output.filename) as f:
  n = len([l for l in f if l.startswith(("ATOM", "HETATM"))
           and l[76:78].strip() == "H"])
assert "mmtbx.ligands.rdkit_utils" not in sys.modules
print("n_h", n)
""" % ([p for p in sys.path if p],)
  result = subprocess.run([sys.executable, "-c", code], capture_output=True,
    universal_newlines=True)
  assert result.returncode == 0, result.stderr[-2000:]
  assert "n_h 48" in result.stdout.splitlines(), result.stdout[-500:]

def test_ccd_atom_table_not_readable():
  if not has_ccd("FCO", "HC0"):
    print("Skipping test_ccd_atom_table_not_readable: no CCD"); return
  status = reduce_hydrogen.ccd_entry_status
  assert (status("FCO").has_entry, status("FCO").n_h) == (True, 0)
  assert status("HC0").n_h == 2
  assert (status("04").has_entry, status("04").n_h) == (False, None)
  # not counted: unknown, not zero
  for edit in (no_type_column, no_atom_loop):
    with edited_ccd("HC0", edit):
      s = status("HC0")
      assert (s.has_entry, s.n_h) == (True, None), (edit, s.n_h)
  # an atom of unknown element: unknown too
  for mark in "?.":
    with edited_ccd("HC0", unknown_h_type(mark)):
      assert status("HC0").n_h is None, mark
      log, description, printed, hd = run_flag("tst_r2_missing_18", pdb_hc0)
    assert "H in its CCD entry" not in log + description, description
    assert "  HC0 A 580: all H (" in description, description
  # the default run stops as before
  with edited_ccd("HC0", no_type_column):
    text = sorry_text("tst_r2_missing_13", pdb_hc0)
  assert text.splitlines()[0] == FIRST + " HC0", text
  # with the flag: no statement on the H of an entry that was not counted
  with edited_ccd("HC0", no_atom_loop):
    log, description, printed, hd = run_flag("tst_r2_missing_14", pdb_hc0)
  assert "HC0" not in hd, hd
  assert "H in its CCD entry" not in log + description, description
  assert "H of its CCD entry" not in log + description, description
  assert "  HC0 A 580: all H (" in description, description

def test_flag_entry_without_h():
  if not has_ccd("FCO"):
    print("Skipping test_flag_entry_without_h: no CCD"); return
  log, description, printed, hd = run_flag("tst_r2_missing_6", pdb_fco)
  assert "FCO" not in hd and len(hd["CYS"]) > 3, hd
  assert ("  FCO L 499: no H in its CCD entry (the CCD entry could not be "
    "used: Explicit valence for atom # 6 O, 3, is greater than permitted)"
    ) in description.splitlines(), description
  for text in (log, description):
    assert "not assessed (CCD entry could not be used): 1 residue" in text, \
      text
    assert "all H" not in text, text
    assert "o dictionary" not in text, text
  assert "The CCD entry of the residue could not be used" in description
  assert ("  CCD entry could not be used, not assessed: 1 residue: FCO L 499 "
    "(no H in its CCD entry)") in log.splitlines(), log
  # the reader's own line: once
  assert printed.count("Could not read CCD entry FCO") == 1, printed

def test_flag_entry_with_h():
  if not has_ccd("HC0"):
    print("Skipping test_flag_entry_with_h: no CCD"); return
  log, description, printed, hd = run_flag("tst_r2_missing_7", pdb_hc0)
  assert "HC0" not in hd, hd
  line = [l for l in description.splitlines()
          if l.startswith("  HC0 A 580: ")]
  assert len(line) == 1 and line[0].startswith("  HC0 A 580: the 2 H of its "
    "CCD entry (the CCD entry could not be used: Explicit valence"), \
    description
  assert "all H" not in description, description
  assert ("  CCD entry could not be used, not assessed: 1 residue: HC0 A 580 "
    "(the 2 H of its CCD entry)") in log.splitlines(), log
  assert printed.count("Could not read CCD entry HC0") == 1, printed
  # a second job in the same process: the reader's line and reason again
  log, description, printed, hd = run_flag("tst_r2_missing_15", pdb_hc0)
  assert printed.count("Could not read CCD entry HC0") == 1, printed
  assert "the 2 H of its CCD entry (the CCD entry could not be used: " \
    "Explicit valence" in description, description

def test_flag_no_entry_keeps_its_wording():
  # with the Cys of pdb_fco: a model without any H to place stops
  cys = "".join(l + "\n" for l in pdb_fco.splitlines() if " CYS " in l)
  log, description, printed, hd = run_flag("tst_r2_missing_8",
    cryst1 + cys + t19.pdb_04)
  assert ("  04 B 1: all H (no dictionary in the library, the CCD or the "
    "restraint files)") in description.splitlines(), description
  assert "not assessed (no dictionary): 1 residue" in log, log
  assert "could not be used" not in log + description

def test_flag_supplied_h_removed_is_named():
  if not has_ccd("EJ2"):
    print("Skipping test_flag_supplied_h_removed_is_named: no CCD"); return
  log, description, printed, hd = run_flag("tst_r2_missing_9", pdb_ej2)
  assert "EJ2" not in hd, hd
  for text in (log, description):
    lines = text.splitlines()
    i = lines.index("Supplied H/D removed from residues without restraints "
      "(keep_existing_H=True keeps them):")
    assert lines[i + 1] == "  EJ2 L 601: H", text
  # kept: no such line, and the hydride is in the output
  log, description, printed, hd = run_flag("tst_r2_missing_10", pdb_ej2,
    keep_existing_H=True)
  assert hd.get("EJ2") == ["H"], hd
  assert "Supplied H/D removed" not in log + description
  assert "its H/D from the input are kept" in description, description

def test_flag_supplied_h_blank_element():
  # the list and the removal agree on what is an H: no element column
  if not has_ccd("EJ2"):
    print("Skipping test_flag_supplied_h_blank_element: no CCD"); return
  log, description, printed, hd = run_flag("tst_r2_missing_16",
    blank_element(pdb_ej2, "H"))
  with open("tst_r2_missing_16_out.pdb") as f:
    names = [l[12:16].strip() for l in f if l[17:20] == "EJ2"]
  assert len(names) == 8 and "H" not in names, names
  for text in (log, description):
    assert "\n  EJ2 L 601: H\n" in text, text

def test_flag_supplied_h_model_ids():
  if not has_ccd("EJ2"):
    print("Skipping test_flag_supplied_h_model_ids: no CCD"); return
  log, description, printed, hd = run_flag("tst_r2_missing_17",
    two_models(pdb_ej2))
  for text in (log, description):
    assert "\n  EJ2 L 601 (model 7): H\n  EJ2 L 601 (model 42): H\n" in \
      text, text

if __name__ == "__main__":
  test_stop_ccd_entry_not_usable()
  test_stop_no_entry()
  test_stop_optimize_says_unchanged()
  test_stop_long_list_is_capped()
  test_stop_names_model_and_altloc()
  test_stop_unl()
  test_flag_unl_advice()
  test_no_rdkit_needed()
  test_ccd_atom_table_not_readable()
  test_flag_entry_without_h()
  test_flag_entry_with_h()
  test_flag_no_entry_keeps_its_wording()
  test_flag_supplied_h_removed_is_named()
  test_flag_supplied_h_blank_element()
  test_flag_supplied_h_model_ids()
  print("OK")
