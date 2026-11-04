"""reduce2 with qmready=True completes the conformers of residues with partial
alternative conformations (iotbx/pdb/altloc_completion.py) before the
atom-complete check: atoms the conformers of a residue share are borrowed
between them, a truncated conformer still needs capping=True.
"""
import os
from six.moves import cStringIO as StringIO
import iotbx.pdb
from iotbx.cli_parser import run_program
from scitbx.matrix import col
from mmtbx.programs import reduce2
from mmtbx.regression.tst_reduce2_capping import run, check_sorry
from mmtbx.regression.model.tst_model_atom_completeness import tripeptide

HEADING = "Altloc completion (atoms borrowed between conformers of a residue):"
H_NOT_PLACED = "qmready: H not placed (see the H placement report above):"
CRYST1 = "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1\n"

def run_log(model_text, prefix, extra_args, fn_out):
  """run, returning the log"""
  fn_in = prefix + ".pdb"
  with open(fn_in, "w") as f:
    f.write(model_text)
  if os.path.exists(fn_out): os.remove(fn_out)
  logger = StringIO()
  run_program(program_class=reduce2.Program, logger=logger,
    args=[fn_in, "overwrite=True", "output.filename=" + fn_out]
      + list(extra_args))
  return logger.getvalue()

def atom_line(line, altloc, occ, xyz=None):
  if xyz is None: xyz = [float(line[30+8*i:38+8*i]) for i in range(3)]
  return "%s%s%s%8.3f%8.3f%8.3f%6.2f%s" % (
    line[:16], altloc, line[17:30], xyz[0], xyz[1], xyz[2], occ, line[60:])

def xyz_of(line):
  return col([float(line[30+8*i:38+8*i]) for i in range(3)])

def tripeptide_lines():
  return tripeptide.splitlines()

def find(lines, resname, name):
  for l in lines:
    if l[17:20] == resname and l[12:16].strip() == name: return l
  raise KeyError((resname, name))

def single_source_model():
  # TYR CZ A/B, OH A only
  lines = []
  for l in tripeptide_lines():
    if l[17:20] == "TYR" and l[12:16] == " CZ ":
      lines.append(atom_line(l, "A", 0.5))
      lines.append(atom_line(l, "B", 0.5))
    elif l[17:20] == "TYR" and l[12:16] == " OH ":
      lines.append(atom_line(l, "A", 1.0))
    else:
      lines.append(l)
  return CRYST1 + "\n".join(lines) + "\n"

def test_single_source_completed():
  model = single_source_model()
  out = "tst_r2ac_1_out.pdb"
  fn_txt = "tst_r2ac_1_out.txt"
  # before altloc completion qmready stopped with "Model must be atom-complete"
  log = run_log(model, "tst_r2ac_1",
    ["qmready=True", "output.description_file_name=" + fn_txt], out)
  h = iotbx.pdb.input(out).construct_hierarchy()
  oh = [a for a in h.atoms() if a.name.strip() == "OH"]
  assert len(oh) == 1 and oh[0].parent().altloc == "", len(oh)
  assert HEADING in log, log
  assert " OH: altloc A -> shared (needed by B)" in log, log
  with open(fn_txt) as f:
    description = f.read()
  assert HEADING in description, description
  assert " OH: altloc A -> shared (needed by B)" in description, description
  # HH on the shared OH, one per conformer of CZ (tst_add_hydrogen_20: the
  # ester rule used to take CZ A and CZ B for two heavy neighbours)
  hh = dict((a.parent().altloc, a) for a in h.atoms()
            if a.name.strip() == "HH")
  assert sorted(hh.keys()) == ["A", "B"], sorted(hh.keys())
  for a in hh.values():
    assert a.distance(oh[0]) < 1.1, a.distance(oh[0])
  for text in (log, description):
    assert H_NOT_PLACED not in text, text

def test_no_heading_without_partial_conformers():
  out = "tst_r2ac_0_out.pdb"
  log = run_log(CRYST1 + tripeptide, "tst_r2ac_0", ["qmready=True"], out)
  assert "Altloc completion" not in log, log

# ASN 2 with CB under A (0.81) and B (0.19, shifted 0.7 A), CG/OD1/ND2 under
# A 0.43, B 0.10, C 0.18 and D 0.20 (C and D: chi1 rotations of A)
def chi1_rotated(lines, angle):
  ca = xyz_of(find(lines, "ASN", "CA"))
  cb = xyz_of(find(lines, "ASN", "CB"))
  axis = (cb - ca).normalize()
  result = {}
  for name in ("CG", "OD1", "ND2"):
    p = xyz_of(find(lines, "ASN", name)) - cb
    result[name] = (p.rotate_around_origin(axis, angle, deg=True) + cb).elems
  return result

def dist(a, b):
  return (col(a) - col(b)).length()

def multi_source_model(occ_d=0.20):
  lines = tripeptide_lines()
  ca = xyz_of(find(lines, "ASN", "CA"))
  cb = xyz_of(find(lines, "ASN", "CB"))
  cg = xyz_of(find(lines, "ASN", "CG"))
  axis = (cb - ca).normalize()
  # CB-B 0.7 A from CB-A, perpendicular to the CA-CB axis, 60 degrees from
  # the direction of CG-A, so that -120 is near it and +120 far
  perp = (cg - cb) - axis * (cg - cb).dot(axis)
  shift = perp.normalize().rotate_around_origin(axis, -60, deg=True) * 0.7
  cb_b = (cb + shift).elems
  c_sites = chi1_rotated(lines, 120)
  # the D angle: -120 unless that fails the bond precondition
  for angle in ([-120] + list(range(-110, -180, -10))
                + list(range(-130, -240, -10))):
    d_sites = chi1_rotated(lines, angle)
    if (dist(d_sites["CG"], cb.elems) < 1.92
        and dist(d_sites["CG"], cb_b) < 1.92): break
  else:
    raise AssertionError("no D angle")
  # both candidates bonded to CG-D, only CB-A to CG-C
  assert dist(c_sites["CG"], cb.elems) < 1.92
  assert dist(c_sites["CG"], cb_b) > 1.92, dist(c_sites["CG"], cb_b)
  assert dist(d_sites["CG"], cb.elems) < 1.92
  assert dist(d_sites["CG"], cb_b) < 1.92
  assert dist(cb_b, ca.elems) < 1.92
  out = []
  for l in lines:
    if l[17:20] != "ASN" or l[12:16].strip() not in ("CB", "CG", "OD1", "ND2"):
      out.append(l)
    elif l[12:16] == " CB ":
      out.append(atom_line(l, "A", 0.81))
      out.append(atom_line(l, "B", 0.19, cb_b))
    else:
      name = l[12:16].strip()
      out.append(atom_line(l, "A", 0.43))
      out.append(atom_line(l, "B", 0.10))
      out.append(atom_line(l, "C", 0.18, c_sites[name]))
      out.append(atom_line(l, "D", occ_d, d_sites[name]))
  return CRYST1 + "\n".join(out) + "\n"

def test_multi_source_copies():
  model = multi_source_model()
  out = "tst_r2ac_2_out.pdb"
  run_log(model, "tst_r2ac_2", ["qmready=True"], out)
  h = iotbx.pdb.input(out).construct_hierarchy()
  asn = [rg for rg in h.residue_groups() if rg.resseq.strip() == "2"][0]
  occ = {}
  hb = {}
  for ag in asn.atom_groups():
    for a in ag.atoms():
      if a.name.strip() == "CB": occ[ag.altloc] = a.occ
      if a.name.strip() in ("HB2", "HB3"):
        hb.setdefault(ag.altloc, []).append(a.name.strip())
  expected = {"A": 0.43, "B": 0.19, "C": 0.18, "D": 0.20}
  assert sorted(occ) == sorted(expected), occ
  for k, v in expected.items():
    assert abs(occ[k] - v) < 1.e-6, (k, occ[k], v)
  for k in expected:
    assert sorted(hb.get(k, [])) == ["HB2", "HB3"], (k, hb)
  # the backbone stays blank
  for ag in asn.atom_groups():
    if ag.altloc == "":
      assert sorted(a.name.strip() for a in ag.atoms()
        if a.element.strip() != "H") == ["C", "CA", "N", "O"]

def test_ambiguous_sorry():
  check_sorry(multi_source_model(0.45), "tst_r2ac_3", ["qmready=True"],
    "Cannot complete conformer D", "ASN", "CB", "altloc A", "altloc B")

def truncated_model():
  lines = tripeptide_lines()
  c1 = xyz_of(find(lines, "GLN", "C"))
  n3 = xyz_of(find(lines, "TYR", "N"))
  axis = (n3 - c1).normalize()
  def rot90(l):
    p = xyz_of(l) - c1
    return (p.rotate_around_origin(axis, 90, deg=True) + c1).elems
  cb = xyz_of(find(lines, "ASN", "CB"))
  ca = xyz_of(find(lines, "ASN", "CA"))
  ca_b = col(rot90(find(lines, "ASN", "CA")))
  assert (ca_b - cb).length() > 1.92, (ca_b - cb).length()
  c_sites = chi1_rotated(lines, 120)
  out = []
  for l in lines:
    name = l[12:16].strip()
    if l[17:20] != "ASN":
      out.append(l)
    elif name in ("N", "CA", "C", "O"):
      out.append(atom_line(l, "A", 0.82))
      out.append(atom_line(l, "B", 0.18, rot90(l)))
    elif name == "CB":
      out.append(atom_line(l, "A", 0.82))
    else:
      out.append(atom_line(l, "A", 0.41))
      out.append(atom_line(l, "C", 0.41, c_sites[name]))
  return CRYST1 + "\n".join(out) + "\n"

def test_truncated_still_needs_capping():
  model = truncated_model()
  check_sorry(model, "tst_r2ac_4", ["qmready=True"],
    "Model must be atom-complete")
  out = "tst_r2ac_4b_out.pdb"
  run(model, "tst_r2ac_4b", ["qmready=True", "capping=True"], fn_out=out)
  h = iotbx.pdb.input(out).construct_hierarchy()
  asn = [rg for rg in h.residue_groups() if rg.resseq.strip() == "2"][0]
  names_b = [a.name.strip() for ag in asn.atom_groups() if ag.altloc == "B"
    for a in ag.atoms()]
  assert "HCB" in names_b and "CB" not in names_b, names_b
  # round trip: the caps of the truncated conformer are not a change
  from iotbx.pdb import altloc_completion
  r = altloc_completion.run(h, check_only=True)
  assert not r.changed, (r.shared, r.copied)
  with open(out) as f:
    capped = f.read()
  log = run_log(capped, "tst_r2ac_4c", ["qmready=True", "capping=True"],
    "tst_r2ac_4c_out.pdb")
  # the truncation of B is reported again (CB, and HCB in A and C), but
  # nothing is shared or copied
  assert HEADING in log, log
  assert " -> shared " not in log and ": copied " not in log, log

def test_constants_agree():
  from mmtbx.hydrogens import capping
  from iotbx.pdb import altloc_completion
  assert altloc_completion.BOND_TOLERANCE == capping.GUESS_BOND_TOLERANCE

def run_all():
  test_constants_agree()
  test_no_heading_without_partial_conformers()
  test_single_source_completed()
  test_multi_source_copies()
  test_ambiguous_sorry()
  test_truncated_still_needs_capping()

if __name__ == "__main__":
  run_all()
  print("OK")
