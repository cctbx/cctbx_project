from __future__ import absolute_import, division, print_function
import os, time
import iotbx.pdb
import iotbx.pdb.hierarchy
import mmtbx.model
from mmtbx.building import split_into_conformers as sic
from libtbx.utils import Sorry, null_out
from libtbx.test_utils import approx_equal
from scitbx import matrix

# Ideal alpha-helix ASKLSEAS, geometry minimized, in a P1 box.
base_str = """\
CRYST1   22.306   22.304   24.857  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.489  14.056  18.803  1.00 10.00           N
ATOM      2  CA  ALA A   1      11.356  15.120  17.815  1.00 10.00           C
ATOM      3  C   ALA A   1      12.115  14.774  16.537  1.00 10.00           C
ATOM      4  O   ALA A   1      11.784  15.262  15.456  1.00 10.00           O
ATOM      5  CB  ALA A   1      11.851  16.440  18.386  1.00 10.00           C
ATOM      6  N   SER A   2      13.137  13.927  16.670  1.00 10.00           N
ATOM      7  CA  SER A   2      13.924  13.535  15.505  1.00 10.00           C
ATOM      8  C   SER A   2      13.263  12.389  14.750  1.00 10.00           C
ATOM      9  O   SER A   2      13.449  12.249  13.535  1.00 10.00           O
ATOM     10  CB  SER A   2      15.340  13.151  15.933  1.00 10.00           C
ATOM     11  OG  SER A   2      16.124  12.767  14.817  1.00 10.00           O
ATOM     12  N   LYS A   3      12.487  11.559  15.450  1.00 10.00           N
ATOM     13  CA  LYS A   3      11.822  10.440  14.790  1.00 10.00           C
ATOM     14  C   LYS A   3      10.506  10.873  14.157  1.00 10.00           C
ATOM     15  O   LYS A   3       9.971  10.178  13.286  1.00 10.00           O
ATOM     16  CB  LYS A   3      11.592   9.305  15.787  1.00 10.00           C
ATOM     17  CG  LYS A   3      12.867   8.755  16.407  1.00 10.00           C
ATOM     18  CD  LYS A   3      12.565   7.632  17.385  1.00 10.00           C
ATOM     19  CE  LYS A   3      13.840   7.082  18.005  1.00 10.00           C
ATOM     20  NZ  LYS A   3      13.560   5.981  18.967  1.00 10.00           N
ATOM     21  N   LEU A   4       9.967  12.018  14.582  1.00 10.00           N
ATOM     22  CA  LEU A   4       8.707  12.492  14.019  1.00 10.00           C
ATOM     23  C   LEU A   4       8.938  13.299  12.747  1.00 10.00           C
ATOM     24  O   LEU A   4       8.004  13.522  11.968  1.00 10.00           O
ATOM     25  CB  LEU A   4       7.950  13.325  15.054  1.00 10.00           C
ATOM     26  CG  LEU A   4       7.577  12.612  16.355  1.00 10.00           C
ATOM     27  CD1 LEU A   4       6.833  13.554  17.290  1.00 10.00           C
ATOM     28  CD2 LEU A   4       6.748  11.369  16.070  1.00 10.00           C
ATOM     29  N   SER A   5      10.175  13.746  12.519  1.00 10.00           N
ATOM     30  CA  SER A   5      10.467  14.525  11.320  1.00 10.00           C
ATOM     31  C   SER A   5      10.845  13.623  10.152  1.00 10.00           C
ATOM     32  O   SER A   5      10.641  13.987   8.988  1.00 10.00           O
ATOM     33  CB  SER A   5      11.582  15.530  11.606  1.00 10.00           C
ATOM     34  OG  SER A   5      11.883  16.299  10.455  1.00 10.00           O
ATOM     35  N   GLU A   6      11.396  12.442  10.441  1.00 10.00           N
ATOM     36  CA  GLU A   6      11.783  11.528   9.371  1.00 10.00           C
ATOM     37  C   GLU A   6      10.582  10.750   8.847  1.00 10.00           C
ATOM     38  O   GLU A   6      10.592  10.279   7.704  1.00 10.00           O
ATOM     39  CB  GLU A   6      12.869  10.573   9.863  1.00 10.00           C
ATOM     40  CG  GLU A   6      14.140  11.265  10.329  1.00 10.00           C
ATOM     41  CD  GLU A   6      15.193  10.287  10.810  1.00 10.00           C
ATOM     42  OE1 GLU A   6      14.931   9.066  10.780  1.00 10.00           O
ATOM     43  OE2 GLU A   6      16.284  10.738  11.218  1.00 10.00           O
ATOM     44  N   ALA A   7       9.539  10.605   9.667  1.00 10.00           N
ATOM     45  CA  ALA A   7       8.356   9.871   9.229  1.00 10.00           C
ATOM     46  C   ALA A   7       7.349  10.795   8.555  1.00 10.00           C
ATOM     47  O   ALA A   7       6.392  10.326   7.926  1.00 10.00           O
ATOM     48  CB  ALA A   7       7.715   9.152  10.416  1.00 10.00           C
ATOM     49  N   SER A   8       7.546  12.104   8.672  1.00 10.00           N
ATOM     50  CA  SER A   8       6.641  13.072   8.064  1.00 10.00           C
ATOM     51  C   SER A   8       7.150  13.514   6.696  1.00 10.00           C
ATOM     52  O   SER A   8       6.415  13.481   5.709  1.00 10.00           O
ATOM     53  CB  SER A   8       6.462  14.287   8.978  1.00 10.00           C
ATOM     54  OG  SER A   8       5.584  15.236   8.398  1.00 10.00           O
TER
END
"""

# ------------------------------------------------------------------ helpers

def get_hierarchy(pdb_str=base_str):
  return iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()

def get_cs(pdb_str=base_str):
  return iotbx.pdb.input(source_info=None, lines=pdb_str).crystal_symmetry()

def get_model(h, cs=None):
  """Unprocessed model from a hierarchy (base cell unless cs is given)."""
  if(cs is None): cs = get_cs()
  h.atoms().reset_serial()
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=h,
    crystal_symmetry=cs, log=null_out())

def fcalc(model, d_min=1.5):
  return model.get_xray_structure().structure_factors(
    d_min=d_min, algorithm="direct").f_calc()

def r_factor(f1, f2):
  a1, a2 = abs(f1).data(), abs(f2).data()
  assert a1.size() == a2.size()
  from scitbx.array_family import flex
  return flex.sum(flex.abs(a1-a2))/flex.sum(a1)

def per_conformer_counts(h, polymer_only=False):
  result = {}
  for rg in h.residue_groups():
    if(polymer_only and not sic._is_polymer(rg)): continue
    for ag in rg.atom_groups():
      result[ag.altloc] = result.get(ag.altloc, 0) + ag.atoms_size()
  return result

def occupancy_sums(h):
  result = {}
  for a in h.atoms():
    ag = a.parent()
    rg = ag.parent()
    key = (rg.parent().id, rg.resseq_as_int(), rg.icode, a.name)
    result[key] = result.get(key, 0) + a.occ
  return result

def occ_of(h, resseq, name):
  """dict altloc -> occupancy of atom name in residue resseq"""
  result = {}
  for a in h.atoms():
    ag = a.parent()
    if(ag.parent().resseq_as_int() == resseq and a.name.strip() == name):
      assert ag.altloc not in result
      result[ag.altloc] = a.occ
  return result

def xyz_of(h, resseq, name):
  result = {}
  for a in h.atoms():
    ag = a.parent()
    if(ag.parent().resseq_as_int() == resseq and a.name.strip() == name):
      result[ag.altloc] = a.xyz
  return result

def residue_group(h, resseq, chain_id="A"):
  for chain in h.only_model().chains():
    if(chain.id != chain_id): continue
    for rg in chain.residue_groups():
      if(rg.resseq_as_int() == resseq): return rg
  raise RuntimeError("no residue %d" % resseq)

def shift(dx, dy, dz):
  return lambda xyz: (xyz[0]+dx, xyz[1]+dy, xyz[2]+dz)

def same():
  return lambda xyz: xyz

def rotation(h, resseq, name_1, name_2, angle):
  """Rotation by angle (deg) about the axis through two atoms of a residue."""
  rg = residue_group(h, resseq)
  p = dict([(a.name.strip(), a.xyz) for a in rg.atoms()])
  p1, p2 = p[name_1], p[name_2]
  return lambda xyz: matrix.rotate_point_around_axis(
    axis_point_1=p1, axis_point_2=p2, point=xyz, angle=angle, deg=True)

def split_atoms(h, resseq, names, states, resnames=None):
  """
  Give atoms of a residue alternative conformations. names: atom names or
  None for all atoms. states: list of (label, occupancy, move) with move a
  function xyz -> xyz. The atoms must be blank.
  """
  rg = residue_group(h, resseq)
  ag = [g for g in rg.atom_groups() if g.altloc.strip() == ""][0]
  if(names is None): names = [a.name.strip() for a in ag.atoms()]
  atoms = [a for a in ag.atoms() if a.name.strip() in names]
  assert len(atoms) == len(names)
  for label, occ, move in states:
    new = iotbx.pdb.hierarchy.atom_group(altloc=label, resname=ag.resname)
    for a in atoms:
      c = a.detached_copy()
      c.set_xyz(move(a.xyz))
      c.set_occ(occ)
      new.append_atom(c)
    rg.append_atom_group(new)
  for a in atoms: ag.remove_atom(a)
  if(ag.atoms_size() == 0): rg.remove_atom_group(ag)

def ser_split(h, resseq, occ_a, occ_b, angle=120):
  """SER OG in two positions: A as is, B rotated about CA-CB."""
  split_atoms(h, resseq, ["OG"], [
    ("A", occ_a, same()),
    ("B", occ_b, rotation(h, resseq, "CA", "CB", angle))])

def run_sic(h, cs=None, **kw):
  return sic.run(model=get_model(h, cs=cs), **kw)

def check_common(result, n_conformers):
  h = result.hierarchy
  counts = per_conformer_counts(h, polymer_only=True)
  assert len(counts) == n_conformers, counts
  assert "" not in counts and " " not in counts, counts
  assert len(set(counts.values())) == 1, counts
  return counts

def expect_sorry(func, text=None):
  try: func()
  except Sorry as e:
    if(text is not None): assert text in str(e), str(e)
    return str(e)
  raise AssertionError("Sorry expected")

# ------------------------------------------------------------------- task 1

def exercise_slot_states():
  assert sic.slot_states([0.3,0.7], 4) == [0,0,1,1]
  assert sic.slot_states([0.3,0.7], 3) == [0,1,1]
  assert sic.slot_states([0.7,0.3], 2) == [1,0]
  assert sic.slot_states([1.0], 3) == [0,0,0]
  assert sic.slot_states([0.2,0.3,0.5], 4) == [0,1,2,2]

def exercise_normalize():
  assert approx_equal(sic.normalized_state_occupancies([0.5,0.4]), [0.55,0.45])
  assert approx_equal(sic.normalized_state_occupancies([0.5,0.4], fill=False),
    [0.5,0.4])
  assert approx_equal(sic.normalized_state_occupancies([0.6,0.6], fill=False),
    [0.5,0.5])
  assert approx_equal(sic.normalized_state_occupancies([0.6,0.6]), [0.5,0.5])
  assert approx_equal(sic.normalized_state_occupancies([0.4,0.6]), [0.4,0.6])
  assert approx_equal(sic.normalized_state_occupancies([0.6]), [0.6])
  assert approx_equal(sic.normalized_state_occupancies([0.,0.]), [0.,0.])

hoh_str = """\
HETATM   55  O   HOH A 101       4.000   4.000   4.000  1.00 20.00           O
HETATM   56  O   HOH A 102       4.000   8.000   4.000  1.00 20.00           O
ATOM     57  H   ALA A   1      11.000  13.200  18.600  1.00 10.00           H
ATOM     58  HA  ALA A   1      10.300  15.200  17.600  1.00 10.00           H
END
"""

def exercise_strip():
  h = get_hierarchy(base_str.replace("TER\nEND\n", hoh_str))
  ser_split(h, 2, 0.4, 0.6)
  result = run_sic(h)
  rh = result.hierarchy
  assert rh.atom_selection_cache().selection(
    "element H or element D or water").count(True) == 0
  assert result.report.n_h_removed == 2
  assert result.report.n_water_removed == 2
  check_common(result, 2)

def exercise_no_altlocs():
  h = get_hierarchy(base_str.replace("TER\nEND\n", hoh_str))
  model = get_model(h)
  result = sic.run(model=model)
  assert result.model is model
  assert result.report.nothing_done
  assert result.hierarchy.atoms_size() == 58

def exercise_equal_counts_and_fcalc():
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  counts = check_common(result, 2)
  assert counts == {"A":54, "B":54}, counts
  rh = result.hierarchy
  assert approx_equal(sorted(occ_of(rh, 2, "OG").values()), [0.4,0.6])
  assert approx_equal(occ_of(rh, 2, "CA"), {"A":0.5, "B":0.5})
  assert approx_equal(occ_of(rh, 6, "CD"), {"A":0.5, "B":0.5})
  for v in occupancy_sums(rh).values(): assert approx_equal(v, 1)
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  assert result.report.labels == ["A","B"]
  assert not result.report.nothing_done

def exercise_four_conformers():
  h = get_hierarchy()
  split_atoms(h, 2, ["OG"], [
    ("A", 0.1, same()),
    ("B", 0.2, rotation(h, 2, "CA", "CB", 90)),
    ("C", 0.3, rotation(h, 2, "CA", "CB", 180)),
    ("D", 0.4, rotation(h, 2, "CA", "CB", 270))])
  ser_split(h, 8, 0.3, 0.7)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  check_common(result, 4)
  rh = result.hierarchy
  assert approx_equal(sorted(occ_of(rh, 8, "OG").values()),
    [0.15,0.15,0.35,0.35])
  assert approx_equal(sorted(occ_of(rh, 2, "OG").values()), [0.1,0.2,0.3,0.4])
  assert approx_equal(sorted(occ_of(rh, 4, "CA").values()), [0.25]*4)
  for v in occupancy_sums(rh).values(): assert approx_equal(v, 1)
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_three_conformers():
  h = get_hierarchy()
  split_atoms(h, 2, ["OG"], [
    ("A", 0.2, same()),
    ("B", 0.3, rotation(h, 2, "CA", "CB", 120)),
    ("C", 0.5, rotation(h, 2, "CA", "CB", 240))])
  ser_split(h, 8, 0.3, 0.7)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  check_common(result, 3)
  assert approx_equal(sorted(occ_of(result.hierarchy, 8, "OG").values()),
    [0.3,0.35,0.35])
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_partial():
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  rg = residue_group(h, 6)
  for a in rg.atoms():
    if(a.name.strip() in ["OE1","OE2"]): a.set_occ(0.5)  # blank, partial
    if(a.name.strip() == "CD"): a.set_occ(0)             # blank, zero
  # single label, partial
  split_atoms(h, 5, ["OG"], [("A", 0.6, same())])
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  check_common(result, 2)
  rh = result.hierarchy
  assert approx_equal(occ_of(rh, 6, "OE1"), {"A":0.25, "B":0.25})
  assert approx_equal(occ_of(rh, 6, "CD"), {"A":0, "B":0})
  assert approx_equal(occ_of(rh, 5, "OG"), {"A":0.3, "B":0.3})
  x = xyz_of(rh, 5, "OG")
  assert approx_equal(x["A"], x["B"])
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  assert len(result.report.occupancy_changes) == 0

def exercise_sum_not_one():
  def model():
    h = get_hierarchy()
    ser_split(h, 2, 0.5, 0.4)
    ser_split(h, 5, 0.7, 0.6)
    return h
  # default: sums below 1 stay, sums above 1 are scaled down
  f0 = fcalc(get_model(model()))
  result = run_sic(model())
  rh = result.hierarchy
  assert approx_equal(sorted(occ_of(rh, 2, "OG").values()), [0.4,0.5])
  assert approx_equal(sorted(occ_of(rh, 5, "OG").values()),
    [0.6/1.3, 0.7/1.3])
  changes = result.report.occupancy_changes
  assert len(changes) == 1, changes
  assert "OG" in changes[0][0] and "SER A   5" in changes[0][0], changes
  assert approx_equal(changes[0][1:], (1.3,1.0))
  assert r_factor(f0, fcalc(result.model)) > 1.e-4
  # without the atom above 1 nothing changes
  h = get_hierarchy()
  ser_split(h, 2, 0.5, 0.4)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  assert result.report.occupancy_changes == []
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  # on request the missing portion is added equally to each state
  result = run_sic(model(), fill_occupancy_to_one=True)
  rh = result.hierarchy
  assert approx_equal(sorted(occ_of(rh, 2, "OG").values()), [0.45,0.55])
  assert approx_equal(sorted(occ_of(rh, 5, "OG").values()),
    [0.6/1.3, 0.7/1.3])
  changes = result.report.occupancy_changes
  assert len(changes) == 2, changes
  assert "OG" in changes[0][0] and "SER" in changes[0][0], changes
  assert approx_equal([c[1:] for c in changes], [(0.9,1.0),(1.3,1.0)])

def exercise_microheterogeneity():
  # residue 5: SER in A, ALA (no OG) in B
  h = get_hierarchy()
  rg = residue_group(h, 5)
  ag = rg.only_atom_group()
  ag.altloc = "A"
  for a in ag.atoms(): a.set_occ(0.6)
  ala = iotbx.pdb.hierarchy.atom_group(altloc="B", resname="ALA")
  for a in ag.atoms():
    if(a.name.strip() == "OG"): continue
    c = a.detached_copy()
    c.set_xyz((a.xyz[0]+0.3, a.xyz[1], a.xyz[2]))
    c.set_occ(0.4)
    ala.append_atom(c)
  rg.append_atom_group(ala)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  counts = per_conformer_counts(result.hierarchy)
  assert sorted(counts.values()) == [53,54], counts
  names = {}
  for g in residue_group(result.hierarchy, 5).atom_groups():
    names[g.resname] = sorted([a.name.strip() for a in g.atoms()])
    assert len(set([round(a.occ,3) for a in g.atoms()])) == 1
  assert names == {"SER":["C","CA","CB","N","O","OG"],
                   "ALA":["C","CA","CB","N","O"]}, names
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_microheterogeneity_blank_atoms():
  # as above, but N and C are blank (under the name SER): both states get them
  h = get_hierarchy()
  rg = residue_group(h, 5)
  ag = rg.only_atom_group()
  ser = iotbx.pdb.hierarchy.atom_group(altloc="A", resname="SER")
  ala = iotbx.pdb.hierarchy.atom_group(altloc="B", resname="ALA")
  for a in ag.atoms():
    if(a.name.strip() in ["N","C"]): continue
    for g, occ, dx in [(ser, 0.6, 0), (ala, 0.4, 0.3)]:
      if(g is ala and a.name.strip() == "OG"): continue
      c = a.detached_copy()
      c.set_xyz((a.xyz[0]+dx, a.xyz[1], a.xyz[2]))
      c.set_occ(occ)
      g.append_atom(c)
    ag.remove_atom(a)
  rg.append_atom_group(ser)
  rg.append_atom_group(ala)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  assert result.report.bridge_atoms == []
  assert result.report.broken_bonds == []
  counts = per_conformer_counts(result.hierarchy)
  assert sorted(counts.values()) == [53,54], counts
  for g in residue_group(result.hierarchy, 5).atom_groups():
    names = sorted([a.name.strip() for a in g.atoms()])
    assert "N" in names and "C" in names, names
    assert len(names) == {"SER":6, "ALA":5}[g.resname]
  assert approx_equal(sorted(occ_of(result.hierarchy, 5, "N").values()),
    [0.5,0.5])
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_microheterogeneity_two_residue_groups():
  # 1EJG PRO/SER 22 pattern: a blank atom under one residue name, then the
  # two residues. Reading the file gives two residue groups of one number.
  lines = []
  for l in base_str.splitlines():
    if(l[17:26] != "SER A   5"):
      lines.append(l)
      continue
    if(l[12:16] == " N  "):
      lines.append(l)
    else:
      lines.append(l[:16]+"A"+l[17:54]+"  0.60"+l[60:])
  i = max([i for i, l in enumerate(lines) if l[17:26] == "SER A   5"])
  ala = []
  for l in base_str.splitlines():
    if(l[17:26] == "SER A   5" and l[12:16] not in [" N  ", " OG "]):
      x = float(l[30:38]) + 0.3
      ala.append(l[:16]+"BALA"+l[20:30]+"%8.3f"%x+l[38:54]+"  0.40"+l[60:])
  lines = lines[:i+1] + ala + lines[i+1:]
  h = get_hierarchy("\n".join(lines))
  assert len([rg for rg in h.residue_groups() if rg.resseq_as_int()==5]) == 2
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  rh = result.hierarchy
  assert len([rg for rg in rh.residue_groups() if rg.resseq_as_int()==5]) == 1
  assert result.report.bridge_atoms == []
  assert result.report.broken_bonds == []
  assert result.report.clashes == []
  counts = per_conformer_counts(rh)
  assert sorted(counts.values()) == [53,54], counts
  groups = residue_group(rh, 5).atom_groups()
  # a position shared by two residue types keeps the labels of the input
  assert sorted([(g.altloc, g.resname) for g in groups]) == [
    ("A","SER"), ("B","ALA")], [(g.altloc, g.resname) for g in groups]
  for g in groups:
    assert "N" in [a.name.strip() for a in g.atoms()]
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_fill_keeps_total():
  # SER 2: CB in A and B (B too far from OG to share it), OG in A only at
  # occupancy 0.4. State B gets a copy of OG, the total stays 0.4.
  h = get_hierarchy()
  split_atoms(h, 2, ["CB"], [("A", 0.4, same()), ("B", 0.6, shift(-2,0,2))])
  rg = residue_group(h, 2)
  blank, ag_a, ag_b = rg.atom_groups()
  og = [a for a in blank.atoms() if a.name.strip() == "OG"][0]
  c = og.detached_copy()
  c.set_occ(0.4)
  ag_a.append_atom(c)
  blank.remove_atom(og)
  f0 = fcalc(get_model(h.deep_copy()))
  # state B has no OG and is too far from that of A, and from the blank CA:
  # stop, listing both bonds
  text = expect_sorry(lambda: run_sic(h.deep_copy(), clash_distance=0.5),
    "cannot be repaired automatically")
  lines = [l for l in text.splitlines() if l.count("SER A   2") == 2]
  assert len(lines) == 2, text
  assert len([l for l in lines if " OG " in l and " CB " in l]) == 1, text
  assert len([l for l in lines if " CA " in l and " CB " in l]) == 1, text
  # the 1US0 GLU 64 pattern: OG without a label, 2 A from CB of state B
  h2 = get_hierarchy()
  split_atoms(h2, 2, ["CB"], [("A", 0.4, same()), ("B", 0.6, shift(-2,0,2))])
  text = expect_sorry(lambda: run_sic(h2, clash_distance=0.5),
    "cannot be repaired automatically")
  lines = [l for l in text.splitlines() if l.count("SER A   2") == 2]
  assert len(lines) == 2, text
  assert len([l for l in lines if " OG " in l and " CB " in l]) == 1, text
  assert len([l for l in lines if " CA " in l and " CB " in l]) == 1, text
  # with a threshold that lets it pass, the copy takes part of the occupancy
  result = run_sic(h, clash_distance=0.5, bond_stretch_threshold=5.0)
  check_common(result, 2)
  assert len(result.report.filled) == 1, result.report.filled
  assert approx_equal(sum(occ_of(result.hierarchy, 2, "OG").values()), 0.4)
  assert approx_equal(sum(occ_of(result.hierarchy, 2, "CB").values()), 1.0)
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_microheterogeneity_three_states():
  # residue 5: ALA in A (0.4), SER in B and C (0.3 each). OG exists in two
  # of three states: its occupancies must not be raised to add up to 1.
  h = get_hierarchy()
  rg = residue_group(h, 5)
  ag = rg.only_atom_group()
  groups = [iotbx.pdb.hierarchy.atom_group(altloc=l, resname=r)
    for l, r in [("A","ALA"), ("B","SER"), ("C","SER")]]
  for g, occ, dx in zip(groups, [0.4,0.3,0.3], [0,0.3,0.6]):
    for a in ag.atoms():
      if(g.resname == "ALA" and a.name.strip() == "OG"): continue
      c = a.detached_copy()
      c.set_xyz((a.xyz[0]+dx, a.xyz[1], a.xyz[2]))
      c.set_occ(occ)
      g.append_atom(c)
  rg.remove_atom_group(ag)
  for g in groups: rg.append_atom_group(g)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h, clash_distance=0.5)
  assert approx_equal(sorted(occ_of(result.hierarchy, 5, "OG").values()),
    [0.3,0.3])
  assert approx_equal(sum(occ_of(result.hierarchy, 5, "CA").values()), 1.0)
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_more_states_than_labels():
  # one label carrying two residue types, and a second label: three states
  h = get_hierarchy()
  rg = residue_group(h, 5)
  ag = rg.only_atom_group()
  for l, r, dx in [("A","SER",0), ("A","ALA",0.3), ("B","SER",0.6)]:
    g = iotbx.pdb.hierarchy.atom_group(altloc=l, resname=r)
    for a in ag.atoms():
      if(r == "ALA" and a.name.strip() == "OG"): continue
      c = a.detached_copy()
      c.set_xyz((a.xyz[0]+dx, a.xyz[1], a.xyz[2]))
      c.set_occ(0.33)
      g.append_atom(c)
    rg.append_atom_group(g)
  rg.remove_atom_group(ag)
  s = expect_sorry(lambda: run_sic(h), "more states")
  assert "A   5" in s, s

def exercise_ligand_numbered_after_chain():
  # a ligand numbered 9 in chain A with an atom named N, far from residue 8:
  # no bond is expected between them
  lig = "HETATM   55  N   NH4 A   9       4.000   4.000   4.000  1.00 20.00" \
    "           N\nEND\n"
  h = get_hierarchy(base_str.replace("TER\nEND\n", lig))
  ser_split(h, 2, 0.4, 0.6)
  result = run_sic(h)
  assert result.report.broken_bonds == [], result.report.broken_bonds
  assert result.report.bridge_atoms == []

def exercise_multi_model():
  lines = base_str.splitlines()
  atoms = [l for l in lines if l.startswith("ATOM")]
  s = "\n".join([lines[0], "MODEL        1"] + atoms + ["ENDMDL",
    "MODEL        2"] + atoms + ["ENDMDL", "END"])
  h = get_hierarchy(s)
  expect_sorry(lambda: run_sic(h), "One model is required")

def exercise_unknown_residue():
  lig = """\
HETATM   55  C1  XQZ A 101       4.000   4.000   4.000  1.00 20.00           C
HETATM   56  C2  XQZ A 101       5.500   4.000   4.000  1.00 20.00           C
END
"""
  h = get_hierarchy(base_str.replace("TER\nEND\n", lig))
  ser_split(h, 2, 0.4, 0.6)
  expect_sorry(lambda: run_sic(h), "restraints (CIF) file")

def exercise_ambiguous_completion():
  # 1R6J SER 205 pattern made ambiguous: CB in A, B; OG in A, B, C, D
  h = get_hierarchy()
  rg = residue_group(h, 2)
  p = dict([(a.name.strip(), a.xyz) for a in rg.atoms()])
  split_atoms(h, 2, ["CB"], [("A", 0.81, same()), ("B", 0.19, shift(0,0.7,0))])
  rg = residue_group(h, 2)
  ag = rg.atom_groups()[0]
  og = [a for a in ag.atoms() if a.name.strip() == "OG"][0]
  for label, occ, move in [("A",0.43,same()), ("B",0.10,shift(0,0.7,0)),
      ("C",0.18,shift(0.2,0.35,0)), ("D",0.45,shift(-0.2,0.35,0.1))]:
    g = [g for g in rg.atom_groups() if g.altloc == label]
    if(len(g) == 0):
      g = iotbx.pdb.hierarchy.atom_group(altloc=label, resname="SER")
      rg.append_atom_group(g)
    else: g = g[0]
    c = og.detached_copy()
    c.set_xyz(move(og.xyz))
    c.set_occ(occ)
    g.append_atom(c)
  ag.remove_atom(og)
  expect_sorry(lambda: run_sic(h), "could come from")

def exercise_special_position():
  # P2, ion on the twofold axis along b through the origin
  cs_str = "CRYST1   40.000   22.304   40.000  90.00  90.00  90.00 P 1 2 1"
  ion = """\
HETATM   55 ZN    ZN A 101       0.000   5.000   0.000  0.50 20.00          ZN
END
"""
  s = base_str.replace(base_str.splitlines()[0], cs_str).replace(
    "TER\nEND\n", ion)
  cs = get_cs(s)
  h = get_hierarchy(s)
  ser_split(h, 2, 0.4, 0.6)
  f0 = fcalc(get_model(h.deep_copy(), cs=cs))
  result = run_sic(h, cs=cs)
  check_common(result, 2)
  # the ion is not split
  assert approx_equal(occ_of(result.hierarchy, 101, "ZN"), {"":0.5})
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  # an atom 1 A from its own symmetry copy
  model = get_model(get_hierarchy(s.replace("   0.000   5.000", "   0.500   5.000")),
    cs=cs)
  ids = [a.id_str() for a in model.get_hierarchy().atoms()]
  assert sic.self_symmetry_clashes(model, ids[:-1], 2.1) == []
  found = sic.self_symmetry_clashes(model, ids, 2.1)
  assert len(found) == 1 and "ZN" in found[0][0], found
  assert approx_equal(found[0][1], 1.0, 1.e-3)
  assert sic.self_symmetry_clashes(model, ids, 0.9) == []

# ------------------------------------------------------------------- task 2

def remove_atoms(h, resseq, names):
  rg = residue_group(h, resseq)
  for ag in rg.atom_groups():
    for a in ag.atoms():
      if(a.name.strip() in names): ag.remove_atom(a)

def bond_length(h, resseq, name_1, name_2, altloc):
  p1 = xyz_of(h, resseq, name_1)[altloc]
  p2 = xyz_of(h, resseq, name_2)[altloc]
  return abs(matrix.col(p1)-matrix.col(p2))

def all_xyz(h):
  result = {}
  for a in h.atoms():
    result.setdefault(a.id_str()[:9]+a.id_str()[10:], set()).add(
      tuple([round(x,3) for x in a.xyz]))
  return result

def check_originals_did_not_move(h_in, h_out):
  """Every position of the input is a position of the same atom in output."""
  x_in, x_out = all_xyz(h_in), all_xyz(h_out)
  for key in x_in:
    assert x_in[key].issubset(x_out[key]), key

def exercise_bridge(name, d_ideal):
  h = get_hierarchy()
  ser_split(h, 8, 0.4, 0.6)
  expected = xyz_of(h, 3, name)[""]
  remove_atoms(h, 3, [name])
  h_in = h.deep_copy()
  result = run_sic(h)
  check_common(result, 2)
  rh = result.hierarchy
  assert per_conformer_counts(rh) == {"A":54, "B":54}
  x = xyz_of(rh, 3, name)
  assert approx_equal(x["A"], x["B"])
  assert abs(matrix.col(x["A"])-matrix.col(expected)) < 0.3
  before, after = {"CB":("CA","CG"), "CG":("CB","CD")}[name]
  for other, d in zip([before, after], d_ideal):
    assert abs(bond_length(rh, 3, name, other, "A") - d) < 0.05
  assert approx_equal(occ_of(rh, 3, name), {"A":0.5, "B":0.5})
  check_originals_did_not_move(h_in, rh)
  assert len(result.report.bridge_atoms) == 1, result.report.bridge_atoms
  assert name in result.report.bridge_atoms[0]
  assert "LYS" in result.report.bridge_atoms[0]

def exercise_open_end():
  h = get_hierarchy()
  ser_split(h, 8, 0.4, 0.6)
  remove_atoms(h, 3, ["NZ"])                 # open end
  remove_atoms(h, 6, ["CG","CD","OE1","OE2"]) # truncated at CB
  remove_atoms(h, 4, ["O"])                  # carbonyl O is a leaf
  result = run_sic(h)
  assert result.report.bridge_atoms == []
  assert per_conformer_counts(result.hierarchy) == {"A":48, "B":48}

def exercise_bridge_in_altloc():
  # LYS side chain beyond CB in two conformations, CB absent in both
  h = get_hierarchy()
  # chi1 rotation: CB is on the axis, so it is the same in both states
  move = rotation(h, 3, "CA", "CB", 120)
  expected = {"A":xyz_of(h, 3, "CB")[""]}
  expected["B"] = expected["A"]
  split_atoms(h, 3, ["CG","CD","CE","NZ"], [
    ("A", 0.3, same()), ("B", 0.7, move)])
  remove_atoms(h, 3, ["CB"])
  h_in = h.deep_copy()
  result = run_sic(h)
  check_common(result, 2)
  rh = result.hierarchy
  x = xyz_of(rh, 3, "CB")
  occ = occ_of(rh, 3, "CB")
  cg = xyz_of(rh, 3, "CG")
  for label in "AB":
    # which input state this slot holds
    state = [s for s in "AB" if abs(matrix.col(cg[label]) -
      matrix.col(xyz_of(h_in, 3, "CG")[s])) < 1.e-3][0]
    assert abs(matrix.col(x[label])-matrix.col(expected[state])) < 0.3
    assert abs(bond_length(rh, 3, "CB", "CG", label) - 1.52) < 0.05
    assert abs(bond_length(rh, 3, "CB", "CA", label) - 1.53) < 0.05
    assert approx_equal(occ[label], {"A":0.3, "B":0.7}[state])
  check_originals_did_not_move(h_in, rh)
  assert len(result.report.bridge_atoms) == 2

def exercise_bridge_ring():
  # a single HIS residue with ideal geometry, its O in two positions
  from mmtbx.monomer_library import idealized_aa
  def his():
    h = idealized_aa.residue_dict()["his"].deep_copy()
    for a in h.atoms():
      a.set_xyz((a.xyz[0]+10, a.xyz[1]+10, a.xyz[2]+10))
      a.set_b(10)
    rg = h.only_residue_group()
    rg.parent().id = "A"
    for a in rg.atoms():
      if(a.name.strip() == "OXT"): rg.only_atom_group().remove_atom(a)
    resseq = rg.resseq_as_int()
    split_atoms(h, resseq, ["O"], [("A", 0.5, same()),
      ("B", 0.5, shift(0.3,0,0))])
    return h, resseq
  # the ring hangs on CG only: an open end, nothing is added
  h, resseq = his()
  remove_atoms(h, resseq, ["ND1","CD2","CE1","NE2"])
  result = run_sic(h)
  assert result.report.bridge_atoms == [], result.report.bridge_atoms
  # one ring atom missing between two present ones: added
  h, resseq = his()
  expected = xyz_of(h, resseq, "CE1")[""]
  remove_atoms(h, resseq, ["CE1"])
  result = run_sic(h)
  assert len(result.report.bridge_atoms) == 1
  x = xyz_of(result.hierarchy, resseq, "CE1")
  assert abs(matrix.col(x["A"])-matrix.col(expected)) < 0.3

def exercise_missing_residue():
  h = get_hierarchy()
  ser_split(h, 8, 0.4, 0.6)
  rg = residue_group(h, 4)
  rg.parent().remove_residue_group(rg)
  result = run_sic(h)
  assert result.report.bridge_atoms == []
  assert per_conformer_counts(result.hierarchy) == {"A":46, "B":46}

def exercise_bridge_main_chain():
  # N of residue 4 joins residue 3 and CA of residue 4
  h = get_hierarchy()
  ser_split(h, 8, 0.4, 0.6)
  expected = xyz_of(h, 4, "N")[""]
  remove_atoms(h, 4, ["N"])
  result = run_sic(h)
  x = xyz_of(result.hierarchy, 4, "N")
  assert abs(matrix.col(x["A"])-matrix.col(expected)) < 0.3
  assert len(result.report.bridge_atoms) == 1
  # N of the first residue is an open end
  h = get_hierarchy()
  ser_split(h, 8, 0.4, 0.6)
  remove_atoms(h, 1, ["N"])
  assert run_sic(h).report.bridge_atoms == []

# ------------------------------------------------------------------- task 3

def rigid_move(h, scale=1):
  """
  Rotation by 8 deg about an axis through the molecule plus a shift of
  0.8 A times scale.
  """
  p1, p2 = (11.0, 12.0, 6.0), (11.5, 12.5, 19.0)
  def move(xyz):
    x = matrix.rotate_point_around_axis(axis_point_1=p1, axis_point_2=p2,
      point=xyz, angle=8, deg=True)
    return (x[0]+0.5*scale, x[1]+0.5*scale, x[2]+0.38*scale)
  return move

def two_conformers(occ_a=0.4, occ_b=0.6, residues=range(1,9), scale=1):
  """All atoms of the residues in two conformers, B a rigid copy of A."""
  h = get_hierarchy()
  move = rigid_move(h, scale=scale)
  for resseq in residues:
    split_atoms(h, resseq, None, [("A", occ_a, same()), ("B", occ_b, move)])
  return h

def swap_labels(h, resseq, names=None, with_occ=True):
  """
  Exchange position and B (and occupancy if with_occ) of atoms between
  labels A and B.
  """
  rg = residue_group(h, resseq)
  a, b = [dict([(x.name.strip(), x) for x in ag.atoms()])
    for ag in rg.atom_groups() if ag.altloc in ["A", "B"]]
  for name in a:
    if(names is not None and name not in names): continue
    xa, xb = a[name], b[name]
    xyz, occ, bb = xa.xyz, xa.occ, xa.b
    xa.set_xyz(xb.xyz); xa.set_b(xb.b)
    xb.set_xyz(xyz); xb.set_b(bb)
    if(with_occ):
      xa.set_occ(xb.occ)
      xb.set_occ(occ)

def covalent_target(model):
  from cctbx import geometry_restraints
  flags = geometry_restraints.flags.flags(default=True, nonbonded=False)
  return model.get_restraints_manager().geometry.energies_sites(
    sites_cart=model.get_sites_cart(), flags=flags).target

def same_atoms(h1, h2, with_occ=True):
  a1, a2 = h1.atoms(), h2.atoms()
  assert a1.size() == a2.size()
  for x, y in zip(a1, a2):
    assert x.id_str() == y.id_str()
    if(abs(matrix.col(x.xyz)-matrix.col(y.xyz)) > 1.e-3): return False
    if(with_occ and abs(x.occ-y.occ) > 1.e-3): return False
  return True

def exercise_untangle():
  reference = run_sic(two_conformers())
  check_common(reference, 2)
  t0 = covalent_target(reference.model)
  assert sum(reference.report.n_moves.values()) == 0, reference.report.n_moves
  assert reference.report.stretched_bonds == []
  # conformer A is the state of lower occupancy everywhere
  for a in reference.hierarchy.atoms():
    assert approx_equal(a.occ, {"A":0.4, "B":0.6}[a.parent().altloc])
  scrambles = [
    ("residue", [(3, None)]),
    ("tail",    [(5, None), (6, None), (7, None), (8, None)]),
    ("branch",  [(3, ["CG","CD","CE","NZ"])]),
    ("atom",    [(2, ["O"])]),
    ("mixed",   [(1, None), (2, ["OG"]), (4, ["CD1"]), (6, None), (7, None),
                 (8, ["C","O"])])]
  # exchange positions with or without their occupancies
  for with_occ in [True, False]:
    for name, scramble in scrambles:
      if(with_occ and name != "mixed"): continue
      h = two_conformers()
      for resseq, names in scramble:
        swap_labels(h, resseq, names, with_occ=with_occ)
      f0 = fcalc(get_model(h.deep_copy()))
      strained = True
      if(strained and name in ["tail", "mixed"]):
        assert covalent_target(
          run_sic(h.deep_copy(), untangle=False).model) > 10*t0
      result = run_sic(h)
      assert covalent_target(result.model) < t0 + 1.e-3, name
      assert same_atoms(result.hierarchy, reference.hierarchy,
        with_occ=with_occ), name
      assert r_factor(f0, fcalc(result.model)) < 1.e-6
      assert result.report.stretched_bonds == []
      if(strained):
        assert sum(result.report.n_moves.values()) >= 1
        if(name in ["tail", "branch"]):
          assert result.report.n_moves["cut"] >= 1, result.report.n_moves
        for (b0, a0), (b1, a1) in zip(result.report.geometry_start,
                                      result.report.geometry_final):
          assert b1 < b0 and b1 < 0.01, (name, b0, b1)

def exercise_local_target():
  # the change of the target computed from the restraints a move touches
  # equals the change of the target of the whole group of atoms
  h = two_conformers(scale=4)
  for resseq, names in [(1, None), (2, ["OG"]), (4, ["CD1"]), (6, None)]:
    swap_labels(h, resseq, names, with_occ=False)
  model = run_sic(h, untangle=False).model
  a = sic.arrangement(model=model, labels=["A","B"])
  assert len(a.components) == 1
  assert len(a.missing_bonds) > 0
  moves = a._moves(a.components[0])
  kinds = set()
  n = 0
  for kind, members in moves[::5] + moves[-3:]:
    kinds.add(kind)
    t0 = a.component_target(0)
    delta = a.exchange(members, 0, 1)
    t1 = a.component_target(0)
    assert approx_equal(t1-t0, delta, 1.e-6*max(1, abs(t0))), (kind, t1-t0,
      delta)
    if(n % 2 == 0): a.exchange(members, 0, 1) # put back every other one
    n += 1
  assert kinds == set(["cut", "residue", "atom"]), kinds

def exercise_untangle_far_apart():
  # conformers 3.2 A apart: a crossover leaves atoms too far to be linked
  # by the restraints, the bond known from the other conformer is used
  reference = run_sic(two_conformers(scale=4))
  assert reference.report.broken_bonds == []
  for scramble in [[5,6,7,8], [2,3,6]]:
    h = two_conformers(scale=4)
    for resseq in scramble: swap_labels(h, resseq, with_occ=False)
    start = run_sic(h.deep_copy(), untangle=False)
    assert len(start.report.broken_bonds) >= 1
    result = run_sic(h)
    assert same_atoms(result.hierarchy, reference.hierarchy, with_occ=False)
    assert result.report.broken_bonds == []
    assert result.report.stretched_bonds == []

def exercise_stretched_reported():
  # residue 5 has a second state 2.5 A away, its neighbours do not
  h = get_hierarchy()
  split_atoms(h, 5, None, [("A", 0.6, same()), ("B", 0.4, shift(2.5,0,0))])
  result = run_sic(h, clash_distance=0.5)
  # one junction is stretched, the other too long to be a bond at all
  stretched = result.report.stretched_bonds
  assert len(stretched) == 1, stretched
  for id_1, id_2, d, d_ideal in stretched:
    assert "A   5" in id_1 + id_2
    assert d - d_ideal > 0.5
  broken = result.report.broken_bonds
  assert len(broken) == 1, broken
  assert "LEU A   4" in broken[0][0] and "SER A   5" in broken[0][1], broken
  assert broken[0][2] > 2.5
  assert sorted(per_conformer_counts(result.hierarchy).values()) == [54,54]

# ------------------------------------------------------------------- task 4

def state_in_slot(h_in, h_out, resseq, name):
  """dict output altloc -> input altloc of the position held there"""
  x_in, x_out = xyz_of(h_in, resseq, name), xyz_of(h_out, resseq, name)
  result = {}
  for label, xyz in x_out.items():
    match = [l for l in x_in
      if abs(matrix.col(x_in[l])-matrix.col(xyz)) < 1.e-3]
    assert len(match) == 1
    result[label] = match[0]
  return result

def exercise_occupancy_examples():
  # the four examples of the specification: residues 2 and 8, not bonded
  for occ_8, flipped in [((0.4,0.6), False), ((0.6,0.4), True),
                         ((0.2,0.8), False), ((0.8,0.2), True)]:
    h = get_hierarchy()
    ser_split(h, 2, 0.4, 0.6)
    ser_split(h, 8, occ_8[0], occ_8[1])
    h_in = h.deep_copy()
    f0 = fcalc(get_model(h.deep_copy()))
    result = run_sic(h)
    rh = result.hierarchy
    for resseq in [2,8]:
      occ = occ_of(rh, resseq, "OG")
      assert occ["A"] < occ["B"]
    assert state_in_slot(h_in, rh, 2, "OG") == {"A":"A", "B":"B"}
    expected = {"A":"B", "B":"A"} if flipped else {"A":"A", "B":"B"}
    assert state_in_slot(h_in, rh, 8, "OG") == expected
    assert approx_equal(sorted(occ_of(rh, 8, "OG").values()), sorted(occ_8))
    assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_geometry_beats_occupancy():
  # residues 4 and 5 in two rigid states; occupancies disagree in order
  h = get_hierarchy()
  move = shift(0.5, 0.5, 0.4)
  split_atoms(h, 4, None, [("A", 0.40, same()), ("B", 0.60, move)])
  split_atoms(h, 5, None, [("A", 0.55, same()), ("B", 0.45, move)])
  h_in = h.deep_copy()
  result = run_sic(h)
  rh = result.hierarchy
  assert state_in_slot(h_in, rh, 4, "CA") == {"A":"A", "B":"B"}
  assert state_in_slot(h_in, rh, 5, "CA") == {"A":"A", "B":"B"}
  assert approx_equal(occ_of(rh, 4, "C"), {"A":0.40, "B":0.60})
  assert approx_equal(occ_of(rh, 5, "N"), {"A":0.55, "B":0.45})
  for label in "AB":
    d = abs(matrix.col(xyz_of(rh, 4, "C")[label]) -
            matrix.col(xyz_of(rh, 5, "N")[label]))
    assert abs(d-1.33) < 0.05, d
  # same with the labels of residue 5 exchanged in the input
  h = h_in.deep_copy()
  swap_labels(h, 5)
  result = run_sic(h)
  assert same_atoms(result.hierarchy, rh)

def ion_str(xyz, altloc=" ", occ=1.0, serial=55, resseq=101, element="CL"):
  return "HETATM%5d %s  %1s %s A%4d    %8.3f%8.3f%8.3f%6.2f 20.00" \
    "          %s\n" % (serial, element, altloc, element, resseq, xyz[0],
    xyz[1], xyz[2], occ, element)

def ser_and_ion(occ_ser, occ_ion, distance, ion_at_state="A", element="CL"):
  """
  SER 2 OG in states A (as is) and B (rotated); a chloride in two states:
  A at the given distance from OG of the given SER state, B far away.
  Returns the hierarchy.
  """
  h = get_hierarchy()
  ser_split(h, 2, occ_ser[0], occ_ser[1])
  og = matrix.col(xyz_of(h, 2, "OG")[ion_at_state])
  cb = matrix.col(xyz_of(h, 2, "CB")[""])
  near = og + (og-cb).normalize()*distance
  far = (18.0, 4.0, 4.0)
  s = h.as_pdb_string(crystal_symmetry=get_cs()).replace("TER\n","").replace(
    "END\n","")
  if(occ_ion is None):
    s += ion_str(near, occ=1.0, element=element)
  elif(isinstance(occ_ion, tuple) and len(occ_ion) == 1):
    s += ion_str(near, "A", occ_ion[0], 55, element=element)
  elif(isinstance(occ_ion, float)):
    s += ion_str(near, occ=occ_ion, element=element)
  else:
    s += ion_str(near, "A", occ_ion[0], 55, element=element) + \
         ion_str(far, "B", occ_ion[1], 56, element=element)
  return get_hierarchy(s + "END\n")

def exercise_clash_beats_occupancy():
  # SER 2 OG and SER 5 OG in two states each; the low-occupancy states are
  # 1.8 A apart (state A of SER 5 is put there, whatever its geometry)
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  og = matrix.col(xyz_of(h, 2, "OG")["A"])
  cb = matrix.col(xyz_of(h, 2, "CB")[""])
  near = (og + (og-cb).normalize()*1.8).elems
  split_atoms(h, 5, ["OG"], [("A", 0.4, lambda xyz: near), ("B", 0.6, same())])
  h_in = h.deep_copy()
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h, bond_stretch_threshold=10.)
  rh = result.hierarchy
  ser_2 = state_in_slot(h_in, rh, 2, "OG")
  ser_5 = state_in_slot(h_in, rh, 5, "OG")
  assert ser_2 != ser_5, (ser_2, ser_5)
  assert len(result.report.clashes) == 0
  assert result.report.n_clash_permutations >= 1
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  check_common(result, 2)

def exercise_ligand_avoids_clash():
  # low-occupancy states of SER OG and of the ion are 1.8 A apart: the
  # polymer stays in occupancy order and the ion exchanges its labels
  h = ser_and_ion((0.4,0.6), (0.4,0.6), 1.8)
  h_in = h.deep_copy()
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  rh = result.hierarchy
  assert state_in_slot(h_in, rh, 2, "OG") == {"A":"A", "B":"B"}
  assert state_in_slot(h_in, rh, 101, "CL") == {"A":"B", "B":"A"}
  assert approx_equal(occ_of(rh, 101, "CL"), {"A":0.6, "B":0.4})
  assert len(result.report.clashes) == 0
  assert result.report.n_clash_permutations == 0
  assert len(result.report.relabelled) == 2, result.report.relabelled
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  check_common(result, 2)
  # labels are in order in the output
  assert [ag.altloc for ag in residue_group(rh, 101).atom_groups()] == ["A","B"]
  # far apart: the ion stays as it is
  h = ser_and_ion((0.4,0.6), (0.4,0.6), 6.0)
  h_in = h.deep_copy()
  result = run_sic(h)
  assert state_in_slot(h_in, result.hierarchy, 2, "OG") == {"A":"A","B":"B"}
  assert state_in_slot(h_in, result.hierarchy, 101, "CL") == {"A":"A","B":"B"}
  assert result.report.relabelled == []
  # so it does when the polymer is rearranged far from it
  h = ser_and_ion((0.6,0.4), (0.6,0.4), 6.0)
  h_in = h.deep_copy()
  result = run_sic(h)
  assert state_in_slot(h_in, result.hierarchy, 2, "OG") == {"A":"B","B":"A"}
  assert state_in_slot(h_in, result.hierarchy, 101, "CL") == {"A":"A","B":"B"}

def exercise_clash_distance_parameter():
  for clash_distance, together in [(2.2, False), (2.1, True)]:
    h = ser_and_ion((0.4,0.6), (0.4,0.6), 2.15)
    h_in = h.deep_copy()
    result = run_sic(h, clash_distance=clash_distance)
    ser = state_in_slot(h_in, result.hierarchy, 2, "OG")
    ion = state_in_slot(h_in, result.hierarchy, 101, "CL")
    assert (ser == ion) == together, (clash_distance, ser, ion)

def exercise_ligand_follows():
  # state A of the ion is in contact (2.6 A) with state A of SER OG. The
  # polymer is rearranged (A is its low-occupancy state): the ion follows
  h = ser_and_ion((0.6,0.4), (0.6,0.4), 2.6, ion_at_state="A")
  h_in = h.deep_copy()
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  rh = result.hierarchy
  assert state_in_slot(h_in, rh, 2, "OG") == {"A":"B", "B":"A"}
  assert state_in_slot(h_in, rh, 101, "CL") == {"A":"B", "B":"A"}
  assert len(result.report.clashes) == 0
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  # in contact with the other state of SER OG: nothing to follow
  h = ser_and_ion((0.6,0.4), (0.6,0.4), 2.6, ion_at_state="B")
  h_in = h.deep_copy()
  result = run_sic(h)
  assert state_in_slot(h_in, result.hierarchy, 101, "CL") == {"A":"A","B":"B"}

def exercise_ligand_single_copy():
  # the 1R6J pattern: an ion labelled A (0.82) sits where state B (0.18) of
  # SER OG is. One copy of the ion stays; it takes the label of the
  # conformer that holds the other state of the residue
  h = ser_and_ion((0.82,0.18), (0.82,), 1.2, ion_at_state="B")
  h_in = h.deep_copy()
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  rh = result.hierarchy
  assert state_in_slot(h_in, rh, 2, "OG") == {"A":"B", "B":"A"}
  assert approx_equal(occ_of(rh, 2, "OG"), {"A":0.18, "B":0.82})
  assert approx_equal(occ_of(rh, 101, "CL"), {"B":0.82})
  assert approx_equal(xyz_of(rh, 101, "CL")["B"], xyz_of(h_in, 101, "CL")["A"])
  assert len(result.report.clashes) == 0, result.report.clashes
  assert len(result.report.relabelled) == 1
  assert result.report.n_atoms == [54, 55], result.report.n_atoms
  check_common(result, 2)
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  # the result is stable
  r2 = sic.run(model=get_model(rh.deep_copy()))
  assert same_atoms(rh, r2.hierarchy)
  assert r2.report.relabelled == []
  # the polymer already in order: the ion keeps its label
  h = ser_and_ion((0.18,0.82), (0.18,), 1.2, ion_at_state="B")
  result = run_sic(h)
  assert approx_equal(occ_of(result.hierarchy, 101, "CL"), {"A":0.18})
  assert result.report.relabelled == []

def exercise_ligand_blank():
  # blank ion at full occupancy, far from everything: stays as it is
  h = ser_and_ion((0.4,0.6), None, 6.0)
  result = run_sic(h)
  check_common(result, 2)
  assert approx_equal(occ_of(result.hierarchy, 101, "CL"), {"":1.0})
  assert result.report.n_atoms == [55, 55]
  # blank ion at occupancy 0.4 overlaps state B of the split SER: alternates
  # encoded by occupancy. The ion gets the label of the other conformer
  h = ser_and_ion((0.6,0.4), 0.4, 1.2, ion_at_state="B")
  h_in = h.deep_copy()
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  rh = result.hierarchy
  ser = state_in_slot(h_in, rh, 2, "OG")
  occ = occ_of(rh, 101, "CL")
  assert list(occ.keys()) == [[l for l in "AB" if ser[l] == "A"][0]], (ser, occ)
  assert approx_equal(list(occ.values()), [0.4])
  assert len(result.report.clashes) == 0
  assert r_factor(f0, fcalc(result.model)) < 1.e-6
  # a full-occupancy blank ion that clashes is reported, and stays blank
  h = ser_and_ion((0.6,0.4), None, 1.2, ion_at_state="B")
  result = run_sic(h)
  assert len(result.report.clashes) == 1, result.report.clashes
  id_1, id_2, d = result.report.clashes[0]
  assert "CL" in id_1+id_2 and "OG" in id_1+id_2
  assert approx_equal(d, 1.2, 1.e-2)
  assert list(occ_of(result.hierarchy, 101, "CL").keys()) == [""]

eoh_str = """\
HETATM   55  C1 AEOH A 101       4.000   4.000   4.000  0.60 20.00           C
HETATM   56  C2 AEOH A 101       5.500   4.000   4.000  0.60 20.00           C
HETATM   57  O  AEOH A 101       6.000   5.300   4.000  0.60 20.00           O
HETATM   58  O  BEOH A 101       8.500   6.000   4.000  0.30 20.00           O
END
"""

def exercise_ligand_truncated_state():
  # the 3NIR pattern: a state of a ligand with one atom only. Ligands are
  # not completed, not copied and their occupancies are not changed
  h = get_hierarchy(base_str.replace("TER\nEND\n", eoh_str))
  ser_split(h, 2, 0.4, 0.6)
  f0 = fcalc(get_model(h.deep_copy()))
  result = run_sic(h)
  check_common(result, 2)
  rg = residue_group(result.hierarchy, 101)
  assert [(ag.altloc, ag.atoms_size()) for ag in rg.atom_groups()] == [
    ("A",3), ("B",1)]
  assert approx_equal(occ_of(result.hierarchy, 101, "O"), {"A":0.6, "B":0.3})
  assert result.report.occupancy_changes == []
  assert result.report.filled == []
  assert r_factor(f0, fcalc(result.model)) < 1.e-6

def exercise_stop():
  # SER 5 OG, blank at occupancy 0.5, overlaps a blank ion: the atom is
  # copied into every conformer and clashes in each
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  for a in residue_group(h, 5).atoms():
    if(a.name.strip() == "OG"): a.set_occ(0.5)
  og = matrix.col(xyz_of(h, 5, "OG")[""])
  cb = matrix.col(xyz_of(h, 5, "CB")[""])
  s = h.as_pdb_string(crystal_symmetry=get_cs()).replace("TER\n","").replace(
    "END\n","") + ion_str(og + (og-cb).normalize()*1.2) + "END\n"
  h = get_hierarchy(s)
  text = expect_sorry(lambda: run_sic(h), "copied into every conformer")
  assert "CL" in text and "OG" in text and "SER A   5" in text, text
  # the same atom at occupancy 0 is a placeholder, not an alternate encoded
  # by occupancy: the clash is reported, no stop
  h0 = get_hierarchy(s)
  for a in residue_group(h0, 5).atoms():
    if(a.name.strip() == "OG"): a.set_occ(0)
  r0 = run_sic(h0)
  assert len(r0.report.clashes) == 1, r0.report.clashes
  assert approx_equal(list(occ_of(r0.hierarchy, 5, "OG").values()), [0, 0])
  # full occupancy, three conformers: reported, no stop; neither when the
  # result is read from a file, where the three copies add up to 0.99
  h = get_hierarchy()
  split_atoms(h, 2, ["OG"], [
    ("A", 0.2, same()),
    ("B", 0.3, rotation(h, 2, "CA", "CB", 120)),
    ("C", 0.5, rotation(h, 2, "CA", "CB", 240))])
  s = h.as_pdb_string(crystal_symmetry=get_cs()).replace("TER\n","").replace(
    "END\n","") + ion_str(og + (og-cb).normalize()*1.2) + "END\n"
  r1 = run_sic(get_hierarchy(s))
  assert len(r1.report.clashes) == 1, r1.report.clashes
  r2 = sic.run(model=get_model(get_hierarchy(r1.model.model_as_pdb())))
  assert len(r2.report.clashes) == 1, r2.report.clashes
  assert approx_equal(sum(occ_of(r2.hierarchy, 5, "OG").values()), 0.99)

def exercise_metal():
  # a metal 2.0 A from an atom is coordination, not a clash; 1.2 A is
  for distance, n_clashes in [(2.0, 0), (1.2, 1)]:
    h = ser_and_ion((0.6,0.4), None, distance, ion_at_state="B", element="ZN")
    result = run_sic(h)
    assert len(result.report.clashes) == n_clashes, result.report.clashes

def exercise_nproc():
  # two side chains in two states each: two independent groups of atoms
  def model(scrambled):
    h = get_hierarchy()
    for resseq, names in [(3, ["CG","CD","CE","NZ"]),
                          (6, ["CG","CD","OE1","OE2"])]:
      split_atoms(h, resseq, names, [("A", 0.4, same()),
        ("B", 0.6, rotation(h, resseq, "CA", "CB", 120))])
    if(scrambled):
      swap_labels(h, 3, ["CE","NZ"], with_occ=False)
      swap_labels(h, 6, ["CD","OE1","OE2"], with_occ=False)
    return h
  reference = run_sic(model(False))
  assert reference.report.n_components == 2
  assert sum(reference.report.n_moves.values()) == 0
  r1 = run_sic(model(True), nproc=1)
  r2 = run_sic(model(True), nproc=2)
  assert same_atoms(r1.hierarchy, r2.hierarchy)
  # each side chain is whole again (which conformer holds which state is
  # open: occupancies of the scrambled atoms were not exchanged)
  assert covalent_target(r2.model) < covalent_target(reference.model) + 1.e-3
  assert r1.report.n_moves == r2.report.n_moves
  assert sum(r2.report.n_moves.values()) >= 2

def exercise_no_symmetry():
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  h.atoms().reset_serial()
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=h,
    log=null_out())
  result = sic.run(model=model)
  check_common(result, 2)

def exercise_idempotent():
  h = two_conformers()
  swap_labels(h, 3)
  r1 = run_sic(h)
  r2 = sic.run(model=get_model(r1.hierarchy.deep_copy()))
  assert same_atoms(r1.hierarchy, r2.hierarchy)
  assert sum(r2.report.n_moves.values()) == 0
  assert r2.report.n_clash_permutations == 0
  # the same through a file, where coordinates are rounded
  for scramble in [[(1, None), (2, ["OG"]), (4, ["CD1"]), (6, None)], []]:
    h = two_conformers()
    for resseq, names in scramble: swap_labels(h, resseq, names, False)
    r1 = run_sic(h)
    s = r1.model.model_as_pdb()
    r2 = sic.run(model=get_model(get_hierarchy(s)))
    assert same_atoms(r1.hierarchy, r2.hierarchy)
    assert sum(r2.report.n_moves.values()) == 0, r2.report.n_moves
    assert r2.report.n_clash_permutations == 0

# ------------------------------------------------------------------- task 5

def exercise_command_line():
  from libtbx import easy_run
  h = get_hierarchy()
  ser_split(h, 2, 0.4, 0.6)
  swap_labels(h, 2, ["OG"], with_occ=False)
  with open("m.pdb", "w") as f:
    f.write(h.as_pdb_string(crystal_symmetry=get_cs()))
  r = easy_run.fully_buffered(
    "phenix.split_into_conformers m.pdb clash_distance=2.1 "
    "fill_occupancy_to_one=True")
  out = "\n".join(r.stdout_lines)
  assert r.return_code == 0, out + "\n".join(r.stderr_lines)
  assert os.path.isfile("m_split.pdb"), out
  for text in ["Conformers: 2 (A, B)", "Atoms per conformer: A 54, B 54",
               "Clash distance: 2.10", "Bridge atoms added: 0",
               "Occupancy sums changed: 0", "Stretched bonds", "Clashes: 0",
               "m_split.pdb"]:
    assert text in out, (text, out)
  rh = iotbx.pdb.input(file_name="m_split.pdb").construct_hierarchy()
  assert per_conformer_counts(rh) == {"A":54, "B":54}
  occ = occ_of(rh, 2, "OG")
  assert occ["A"] < occ["B"]
  # a model without alternative conformations: nothing written
  with open("n.pdb", "w") as f: f.write(base_str)
  r = easy_run.fully_buffered("phenix.split_into_conformers n.pdb")
  out = "\n".join(r.stdout_lines)
  assert r.return_code == 0
  assert "no alternative conformations" in out, out
  assert not os.path.isfile("n_split.pdb")

def run():
  t0 = time.time()
  folder = "tst_split_into_conformers"
  if(not os.path.isdir(folder)): os.makedirs(folder)
  os.chdir(folder)
  tests = [
    ("exercise_slot_states", exercise_slot_states),
    ("exercise_normalize", exercise_normalize),
    ("exercise_strip", exercise_strip),
    ("exercise_no_altlocs", exercise_no_altlocs),
    ("exercise_equal_counts_and_fcalc", exercise_equal_counts_and_fcalc),
    ("exercise_four_conformers", exercise_four_conformers),
    ("exercise_three_conformers", exercise_three_conformers),
    ("exercise_partial", exercise_partial),
    ("exercise_sum_not_one", exercise_sum_not_one),
    ("exercise_microheterogeneity", exercise_microheterogeneity),
    ("exercise_microheterogeneity_blank_atoms",
      exercise_microheterogeneity_blank_atoms),
    ("exercise_microheterogeneity_two_residue_groups",
      exercise_microheterogeneity_two_residue_groups),
    ("exercise_fill_keeps_total", exercise_fill_keeps_total),
    ("exercise_microheterogeneity_three_states",
      exercise_microheterogeneity_three_states),
    ("exercise_more_states_than_labels", exercise_more_states_than_labels),
    ("exercise_ligand_numbered_after_chain",
      exercise_ligand_numbered_after_chain),
    ("exercise_multi_model", exercise_multi_model),
    ("exercise_unknown_residue", exercise_unknown_residue),
    ("exercise_ambiguous_completion", exercise_ambiguous_completion),
    ("exercise_special_position", exercise_special_position),
    ("exercise_bridge 'CB', (1.53, 1.52)", lambda: exercise_bridge("CB", (1.53, 1.52))),
    ("exercise_bridge 'CG', (1.52, 1.52)", lambda: exercise_bridge("CG", (1.52, 1.52))),
    ("exercise_open_end", exercise_open_end),
    ("exercise_bridge_in_altloc", exercise_bridge_in_altloc),
    ("exercise_bridge_ring", exercise_bridge_ring),
    ("exercise_missing_residue", exercise_missing_residue),
    ("exercise_bridge_main_chain", exercise_bridge_main_chain),
    ("exercise_untangle", exercise_untangle),
    ("exercise_local_target", exercise_local_target),
    ("exercise_untangle_far_apart", exercise_untangle_far_apart),
    ("exercise_stretched_reported", exercise_stretched_reported),
    ("exercise_occupancy_examples", exercise_occupancy_examples),
    ("exercise_geometry_beats_occupancy", exercise_geometry_beats_occupancy),
    ("exercise_clash_beats_occupancy", exercise_clash_beats_occupancy),
    ("exercise_clash_distance_parameter", exercise_clash_distance_parameter),
    ("exercise_ligand_avoids_clash", exercise_ligand_avoids_clash),
    ("exercise_ligand_follows", exercise_ligand_follows),
    ("exercise_ligand_single_copy", exercise_ligand_single_copy),
    ("exercise_ligand_blank", exercise_ligand_blank),
    ("exercise_ligand_truncated_state", exercise_ligand_truncated_state),
    ("exercise_stop", exercise_stop),
    ("exercise_metal", exercise_metal),
    ("exercise_nproc", exercise_nproc),
    ("exercise_no_symmetry", exercise_no_symmetry),
    ("exercise_idempotent", exercise_idempotent),
    ("exercise_command_line", exercise_command_line),]
  for name, test in tests:
    t = time.time()
    test()
    print("  %-36s %5.1f s" % (name, time.time()-t))
  print("Time: %.1f" % (time.time()-t0))
  print("OK")

if(__name__ == "__main__"):
  run()
